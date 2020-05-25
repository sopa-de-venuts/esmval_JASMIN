# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

ESMVAL_ENV_NAME=esmval_primavera

alias sci1='ssh jasmin-sci1'
alias sci2='ssh jasmin-sci2'
alias sci3='ssh jasmin-sci3'
alias xfer='ssh jasmin-xfer1'
alias esmval_env='conda activate /gws/smf/j04/primavera/envs/${ESMVAL_ENV_NAME}'
# User specific aliases and functions

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/users/pcos/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/users/pcos/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/users/pcos/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/users/pcos/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

