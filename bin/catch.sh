
# from https://stackoverflow.com/questions/11027679/capture-stdout-and-stderr-into-different-variables
_CTRL_Z_=$'\cZ'

# SYNTAX:
#     catch_posix STDOUT_VARIABLE STDERR_VARIABLE COMMAND
catch_posix() {
    {
        IFS=$'\n'"${_CTRL_Z_}" read -r -d "${_CTRL_Z_}" "${1}";
        IFS=$'\n'"${_CTRL_Z_}" read -r -d "${_CTRL_Z_}" "${2}";
        (IFS=$'\n'"${_CTRL_Z_}" read -r -d "${_CTRL_Z_}" _ERRNO_; return ${_ERRNO_});
    } <<EOF
$((printf "${_CTRL_Z_}%s${_CTRL_Z_}%d${_CTRL_Z_}" "$(((({ ${3}; echo "${?}" 1>&3-; } | cut -z -d"${_CTRL_Z_}" -f1 | tr -d '\0' 1>&4-) 4>&2- 2>&1- | cut -z -d"${_CTRL_Z_}" -f1 | tr -d '\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
EOF
}
