/*
 * ai_call.c - Setuid wrapper to call HMS AI inference service.
 * Owned by the token holder; other users can run it but never read the token.
 *
 * Compile:  gcc -O2 -Wall -o ai_call ai_call.c -lcurl
 * Install:  chmod 4755 ai_call          (no sudo needed -- you own it)
 *           chmod 600 ~/.smartSlurm/.ai_token
 *
 * Usage:    ai_call "your prompt here"
 *           echo "your prompt" | ai_call
 *
 * Security: The setuid bit makes the binary run with the owner's effective UID,
 *           so it can read the 600-permission token file. Other users cannot.
 *           The token is cleared from memory immediately after use.
 *           The token never appears in process args or environment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <curl/curl.h>

#define AI_BASE_URL    "https://ai-poc.hms.edu/v1/chat/completions"
#define AI_MODEL       "google/gemma-4-31B-it"
#define TOKEN_SUFFIX   "/.smartSlurm/.ai_token"
#define ALLOWED_SUFFIX "/.smartSlurm/allowed_users.txt"
#define MAX_TOKEN      4096
#define MAX_PROMPT     65536

/* Accumulate curl response into a growable buffer */
struct MemBuf {
    char  *data;
    size_t size;
};

static size_t write_cb(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t total = size * nmemb;
    struct MemBuf *buf = (struct MemBuf *)userp;
    char *ptr = realloc(buf->data, buf->size + total + 1);
    if (!ptr) { fprintf(stderr, "Out of memory\n"); return 0; }
    buf->data = ptr;
    memcpy(buf->data + buf->size, contents, total);
    buf->size += total;
    buf->data[buf->size] = '\0';
    return total;
}

/*
 * Check whether the real caller (getuid) is listed in the owner's
 * ~/.smartSlurm/allowed_users.txt (one username per line).
 * The owner is always implicitly allowed.
 * The file is opened as effective UID (owner) since it is chmod 600.
 */
static int check_allowed(void) {
    uid_t caller_uid = getuid();   /* real UID = who actually ran the binary */
    uid_t owner_uid  = geteuid(); /* effective UID = binary owner (setuid)  */

    if (caller_uid == owner_uid) return 0; /* owner is always allowed */

    struct passwd *owner_pw  = getpwuid(owner_uid);
    struct passwd *caller_pw = getpwuid(caller_uid);
    if (!owner_pw) { fprintf(stderr, "Cannot resolve owner home dir\n"); return -1; }

    char path[4096];
    snprintf(path, sizeof(path), "%s%s", owner_pw->pw_dir, ALLOWED_SUFFIX);

    FILE *f = fopen(path, "r"); /* readable because euid == owner */
    if (!f) {
        fprintf(stderr, "Access denied: no allowed_users.txt found.\n");
        fprintf(stderr, "Owner: add usernames to %s\n", path);
        return -1;
    }

    const char *caller_name = caller_pw ? caller_pw->pw_name : NULL;
    char line[256];
    int found = 0;
    while (fgets(line, sizeof(line), f)) {
        /* Strip trailing newline/whitespace */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' '))
            line[--len] = '\0';
        if (len == 0 || line[0] == '#') continue; /* skip blank lines and comments */
        if (caller_name && strcmp(line, caller_name) == 0) { found = 1; break; }
    }
    fclose(f);

    if (!found) {
        fprintf(stderr, "Access denied: '%s' is not in the allowed users list.\n",
                caller_name ? caller_name : "unknown");
        fprintf(stderr, "Owner: add the username to %s\n", path);
        return -1;
    }
    return 0;
}

/*
 * Read token using the effective UID (the binary owner, set by the setuid bit).
 * The effective UID already has permission to open a chmod-600 token file.
 */
static int read_token(char *token, size_t maxlen) {
    uid_t owner_uid = geteuid(); /* effective UID = binary owner (setuid) */
    struct passwd *pw = getpwuid(owner_uid);
    if (!pw) { fprintf(stderr, "Cannot resolve token owner home dir\n"); return -1; }

    char path[4096];
    snprintf(path, sizeof(path), "%s%s", pw->pw_dir, TOKEN_SUFFIX);

    FILE *f = fopen(path, "r"); /* works because euid == file owner */
    if (!f) {
        fprintf(stderr, "Token file not found: %s\n", path);
        fprintf(stderr, "Fix: echo YOUR_TOKEN > %s && chmod 600 %s\n", path, path);
        return -1;
    }

    if (!fgets(token, (int)maxlen, f)) {
        fclose(f); fprintf(stderr, "Token file is empty\n"); return -1;
    }
    fclose(f);

    /* Strip trailing newline/whitespace */
    size_t len = strlen(token);
    while (len > 0 && (token[len-1] == '\n' || token[len-1] == '\r' || token[len-1] == ' '))
        token[--len] = '\0';

    if (len < 10) { fprintf(stderr, "Token looks invalid (too short)\n"); return -1; }
    return 0;
}

/* JSON-escape a string for safe embedding in a JSON string value */
static void json_escape(const char *in, char *out, size_t outlen) {
    size_t j = 0;
    for (size_t i = 0; in[i] && j + 8 < outlen; i++) {
        unsigned char c = (unsigned char)in[i];
        if      (c == '"')  { out[j++] = '\\'; out[j++] = '"';  }
        else if (c == '\\') { out[j++] = '\\'; out[j++] = '\\'; }
        else if (c == '\n') { out[j++] = '\\'; out[j++] = 'n';  }
        else if (c == '\r') { out[j++] = '\\'; out[j++] = 'r';  }
        else if (c == '\t') { out[j++] = '\\'; out[j++] = 't';  }
        else if (c < 0x20)  { j += snprintf(out+j, outlen-j, "\\u%04x", c); }
        else                 { out[j++] = c; }
    }
    out[j] = '\0';
}

/* Print the "content" value from an OpenAI-format JSON response */
static void print_response(const char *json) {
    const char *p = strstr(json, "\"content\":");
    if (!p) { puts(json); return; } /* fallback: dump raw */

    p += strlen("\"content\":");
    while (*p == ' ') p++;
    if (*p != '"') { puts(json); return; }
    p++;

    while (*p && !(*p == '"' && *(p-1) != '\\')) {
        if      (*p == '\\' && *(p+1) == 'n')  { putchar('\n'); p += 2; }
        else if (*p == '\\' && *(p+1) == 't')  { putchar('\t'); p += 2; }
        else if (*p == '\\' && *(p+1) == '"')  { putchar('"');  p += 2; }
        else if (*p == '\\' && *(p+1) == '\\') { putchar('\\'); p += 2; }
        else                                    { putchar(*p);   p++;    }
    }
    putchar('\n');
}

int main(int argc, char *argv[]) {
    if (check_allowed() != 0) return 1;

    char token[MAX_TOKEN];
    if (read_token(token, sizeof(token)) != 0) return 1;

    /* Read prompt from argument(s) or stdin */
    char prompt_raw[MAX_PROMPT];
    if (argc >= 2) {
        prompt_raw[0] = '\0';
        for (int i = 1; i < argc; i++) {
            if (i > 1) strncat(prompt_raw, " ", sizeof(prompt_raw) - strlen(prompt_raw) - 1);
            strncat(prompt_raw, argv[i], sizeof(prompt_raw) - strlen(prompt_raw) - 1);
        }
    } else {
        size_t n = fread(prompt_raw, 1, sizeof(prompt_raw) - 1, stdin);
        prompt_raw[n] = '\0';
    }

    if (strlen(prompt_raw) == 0) {
        fprintf(stderr, "Usage: %s \"prompt\"  OR  echo \"prompt\" | %s\n", argv[0], argv[0]);
        return 1;
    }

    /* JSON-escape the prompt */
    char prompt_esc[MAX_PROMPT * 2];
    json_escape(prompt_raw, prompt_esc, sizeof(prompt_esc));

    /* Build request body */
    char body[MAX_PROMPT * 2 + 256];
    snprintf(body, sizeof(body),
        "{\"model\":\"%s\",\"messages\":[{\"role\":\"user\",\"content\":\"%s\"}]}",
        AI_MODEL, prompt_esc);

    /* Build Authorization header and clear token from memory immediately */
    char auth_hdr[MAX_TOKEN + 32];
    snprintf(auth_hdr, sizeof(auth_hdr), "Authorization: Bearer %s", token);
    memset(token, 0, sizeof(token));

    /* HTTP request via libcurl */
    curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL *curl = curl_easy_init();
    if (!curl) { fprintf(stderr, "curl_easy_init failed\n"); return 1; }

    struct MemBuf resp = {NULL, 0};
    struct curl_slist *hdrs = NULL;
    hdrs = curl_slist_append(hdrs, "Content-Type: application/json");
    hdrs = curl_slist_append(hdrs, auth_hdr);
    memset(auth_hdr, 0, sizeof(auth_hdr)); /* clear header from memory */

    curl_easy_setopt(curl, CURLOPT_URL,           AI_BASE_URL);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER,    hdrs);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS,    body);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_cb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA,     &resp);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT,       120L);

    CURLcode rc = curl_easy_perform(curl);
    curl_slist_free_all(hdrs);
    curl_easy_cleanup(curl);
    curl_global_cleanup();

    if (rc != CURLE_OK) {
        fprintf(stderr, "Request failed: %s\n", curl_easy_strerror(rc));
        free(resp.data);
        return 1;
    }

    if (resp.data) {
        print_response(resp.data);
        free(resp.data);
    }
    return 0;
}
