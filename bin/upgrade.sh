#!/bin/bash

#set -x 

sourcePath=`dirname $0`

cd $sourcePath

TARGET_BRANCH="master"
FILE_TO_KEEP="config/config.txt"

current_branch=$(git symbolic-ref --short HEAD)

echo "Check to see if we are in master branch..."
if [ "$current_branch" != "$TARGET_BRANCH" ]; then
    echo "Currently on branch $current_branch. Please switching to $TARGET_BRANCH first ..."
    exit 1; 
    #git checkout $TARGET_BRANCH
fi

echo "Save local change..."
git stash

echo "Pulling latest changes on $TARGET_BRANCH..."
git pull

echo "Apply stashed changes back"
git stash pop

echo "Checking for conflicts..."
if git diff --name-only --diff-filter=U | grep -q "$FILE_TO_KEEP"; then
    echo "Conflicts detected in $FILE_TO_KEEP. Please resolve them manually."
else
    echo "No conflicts detected. Cleaning up stash."
    git stash drop
fi

echo "Upgrade process completed."