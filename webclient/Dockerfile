FROM webdevops/php-nginx:ubuntu-18.04

WORKDIR /app
COPY . /app


RUN apt update
RUN apt-get install dialog apt-utils
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:inkscape.dev/stable
RUN apt update
RUN apt-get -y install inkscape

RUN apt-get -y install imagemagick && \
  apt-get -y install nodejs && \
  apt-get -y install npm && \
  npm install -g bower && \
  npm install -g grunt-cli && \
  npm install grunt --save-dev && \
  npm install grunt-contrib-clean grunt-contrib-uglify grunt-text-replace grunt-contrib-less grunt-svgmin grunt-contrib-copy grunt-contrib-watch


RUN chmod a+x build.sh
RUN ./build.sh fetch jmol


EXPOSE 80