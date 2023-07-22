#!/usr/bin/env bash

##previous_tag=HEAD
previous_tag=0
FILTER="${1:-.*}"
for current_tag in $(git tag --sort=-creatordate)
do
if [ "$previous_tag" != 0 ];then
    tag_date=$(git log -1 --pretty=format:'%ad' --date=short ${previous_tag})
    printf "##### ${previous_tag} (${tag_date})\n\n"
    git log ${current_tag}...${previous_tag} --pretty=format:'-  %s [view](https://github.com/bigomics/omicsplayground/commit/%H)' --reverse | grep -E -v "(Merge|Style)" | grep -i -E "$FILTER"
    printf "\n\n"
fi
previous_tag=${current_tag}
done
