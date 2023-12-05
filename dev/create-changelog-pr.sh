#!/usr/bin/env bash

previous_tag=0
FILTER="${1:-.*}"
NHEAD="${2:-999}"

for current_tag in $(git tag --sort=-creatordate | head -n ${NHEAD})
do
    if [ "$previous_tag" != 0 ];then
        tag_date=$(git log -1 --pretty=format:'%ad' --date=short ${previous_tag})
        printf "### ${previous_tag} (${tag_date})\n\n"
        for pr in $(git log ${current_tag}...${previous_tag} --pretty="%s" --merges | grep pull | sed -e "s/.*#\([0-9]\+\).*/\1/g" | sort -rn | uniq); do
            PR_title=$(curl --request GET \
            -s \
            --header "Authorization: Bearer ${GITHUB_API_KEY}" \
            --header "X-GitHub-Api-Version: 2022-11-28" \
            --url https://api.github.com/repos/bigomics/omicsplayground/pulls/${pr} | jq -r '.title')

            if [ "$PR_title" != "null" ]; then
                echo "- ${PR_title} [view](https://github.com/bigomics/omicsplayground/pull/${pr})"
            fi
            
           
        done
        printf "\n\n"
    fi
    previous_tag=${current_tag}
done