"""Optional AI support for the SmartSlurm Python tools.

The Python twin of bin/aiHelper.sh. Every function here is safe to call when
AI is unavailable: ai_enabled() returns False and ai_ask() returns None, so
callers keep working exactly as they did before AI existed.

    try:
        import aiHelper
    except ImportError:
        aiHelper = None

The AI token never appears here: it stays inside the setuid ai_call helper.
"""

import os
import subprocess
from pathlib import Path

# ai_call refuses prompts larger than 64K, so stay well under it.
MAX_PROMPT = 48000


def _read_config():
    """Read the SmartSlurm config the same way the shell scripts do."""
    home_config = Path.home() / ".smartSlurm/config/config.txt"
    local_config = Path(__file__).resolve().parent.parent / "config/config.txt"

    if home_config.is_file():
        path = home_config
    elif local_config.is_file():
        path = local_config
    else:
        return {}

    values = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line.startswith("export "):
                    continue
                key, sep, value = line[len("export "):].partition("=")
                if not sep:
                    continue
                value = value.split("#", 1)[0].strip().strip('"').strip("'")
                values[key.strip()] = value.replace("$HOME", str(Path.home()))
    except OSError:
        return {}
    return values


_config = _read_config()


def _setting(name, default=None):
    """Config file first, then the environment.

    The shell scripts source config.txt, which exports over anything already in
    the environment, so the config file wins there too. Keep the same order here.
    """
    value = _config.get(name)
    if not value:
        value = os.environ.get(name)
    return value if value else default


def ai_call_path():
    """Path to the setuid ai_call helper, or None if it cannot be found."""
    explicit = _setting("smartSlurmAiCall")
    if explicit:
        return explicit if os.access(explicit, os.X_OK) else None

    for candidate in (
        Path.home() / ".smartSlurm/bin/ai_call",
        Path(__file__).resolve().parent / "ai_call",
        Path("/home/ld32/.smartSlurm/bin/ai_call"),
    ):
        if os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def ai_enabled():
    """True only when AI is turned on and the helper is really runnable."""
    if _setting("smartSlurmAiEnabled", "yes") == "no":
        return False
    return ai_call_path() is not None


def ai_why_not():
    """One line explaining why AI is unavailable, for the user."""
    if _setting("smartSlurmAiEnabled", "yes") == "no":
        return "AI is disabled by config (smartSlurmAiEnabled=no)."
    if ai_call_path() is None:
        return (
            "AI helper ai_call not found or not executable. Ask your SmartSlurm "
            "admin to install it and add '%s' to the owner's "
            "~/.smartSlurm/allowed_users.txt" % os.environ.get("USER", "your user")
        )
    return (
        "AI call failed. If you see 'Access denied', ask your SmartSlurm admin "
        "to add '%s' to the owner's ~/.smartSlurm/allowed_users.txt"
        % os.environ.get("USER", "your user")
    )


def _timeout(explicit=None):
    if explicit is not None:
        return explicit
    try:
        return int(_setting("smartSlurmAiTimeout", "90"))
    except ValueError:
        return 90


def ai_ask(prompt, timeout=None):
    """Send a prompt to AI and return the answer, or None on any failure.

    Never raises and never blocks longer than the timeout, so a caller can treat
    a None result as simply 'no AI today'.
    """
    if not ai_enabled():
        return None

    helper = ai_call_path()
    if not helper or not prompt.strip():
        return None

    try:
        result = subprocess.run(
            [helper],
            input=prompt[:MAX_PROMPT],
            capture_output=True,
            text=True,
            timeout=_timeout(timeout),
        )
    except (subprocess.TimeoutExpired, OSError):
        return None

    if result.returncode == 0 and result.stdout.strip():
        return result.stdout.strip()
    return None
