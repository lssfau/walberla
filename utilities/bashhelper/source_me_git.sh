###########################################################################################
###################          git prompt                        ############################
###########################################################################################
#
# This script modifies your prompt to better support git repositories
#
# Setup:
#   Insert the following into your .bashrc file:
#   export WALBERLA_SOURCE_DIR=/path/to/walberla/sources
#   source $WALBERLA_SOURCE_DIR/utilities/bashhelper/source_me_git.sh
#
###########################################################################################

WSD=$WALBERLA_SOURCE_DIR

if git --version &> /dev/null && [ -f $WSD/utilities/bashhelper/git-prompt.sh ]; then 
   source $WSD/utilities/bashhelper/git-prompt.sh
   PROMPT_COMMAND='__git_ps1 "\u@\h:\w" "> "'
   GIT_PS1_SHOWDIRTYSTATE=true
   GIT_PS1_SHOWCOLORHINTS=true
fi

