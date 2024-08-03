apt update && apt install -y \
    locales apt-utils software-properties-common \
    libcurl4-gnutls-dev libnode-dev libv8-dev \
    libssl-dev libxml2-dev libjpeg-dev \
    libgl-dev libglu-dev tk-dev libhdf5-dev \
    libgit2-dev libssh2-1-dev libnetcdf-dev \
    libudunits2-dev libgdal-dev libbz2-dev \
    jags emacs git procps htop \
    python3 python3-pip python-is-python3 \
    pdftk vim-tiny less wget gdebi-core \
    pandoc imagemagick libfftw3-dev libglpk-dev \
    libgsl-dev librsvg2-dev libgsl-dev curl \
    libsodium-dev cmake
##libcairo2-dev libxt-dev xvfb \

# Install Chrome
apt-get update && apt-get install -y wget gnupg2
wget -q -O - https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add -
echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main" > /etc/apt/sources.list.d/google-chrome.list
apt-get update && apt-get install -y google-chrome-stable

# Set the locale to UTF8
# sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
# export LANG=en_US.UTF-8
# export LANGUAGE=en_US:en
# export LC_ALL=en_US.UTF-8 

# Install basic stuff and R
apt install -y r-base r-base-dev 

# remove ugly snaps
#snap remove --purge -y firefox
#sudo apt purge snapd
