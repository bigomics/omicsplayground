## run from root

## get version
version=`cat VERSION`
image=bigomics/playground:$version
echo VERSION=$version

## update datasets info
echo '(cd scripts && Rscript update-datasets-info.R)'

##---------------- github ----------------------

echo git pull
echo git tag -a "$version" -m "version $version"
echo git push
echo git push --tags


##---------------- docker ----------------------
echo nohup docker build --no-cache -t $image . > docker.out &
echo docker build -t $image .

## +run in background, remove contained after use, give nice name
echo docker run --rm -d -p 4000:3838 --name=play1 $image

## To save your Docker Image as a tar-archive, you simply type into your terminal:
echo docker save -o ~/playground_$version.tar $image

## publish to Docker Hub
echo docker tag $image bigomics/playground:latest
echo docker push bigomics/playground:$version
echo docker push bigomics/playground:latest
