From ubuntu:18.04

WORKDIR /app
COPY ./webclient /app

RUN apt-get -y install nginx && \
  apt-get -y install php5.6 && \
  apt-get -y install php5.6-cli && \
  apt-get -y install php5.6-fpm && \
  apt-get -y install inkscape && \
  apt-get -y install imagemagick && \
  apt-get -y install php-imagick && \
  apt-get -y install nodejs && \
  apt-get -y install npm && \
  npm install -g bower && \
  npm install -g grunt-cli && \
  npm install grunt --save-dev && \
  npm install grunt-contrib-clean grunt-contrib-uglify grunt-text-replace grunt-contrib-less grunt-svgmin grunt-contrib-copy grunt-contrib-watch

RUN cd /app/webclient
RUN ./build.sh fetch jmol

COPY ./example.com.conf /etc/nginx/conf.d/
RUN nginx -s reload