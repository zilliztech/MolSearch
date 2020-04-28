#docker build -f Dockerfile.os.build -t molsearch:base .
FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
RUN sed -i "s/archive.ubuntu.com/mirrors.aliyun.com/g" /etc/apt/sources.list

RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y nginx php-fpm
RUN echo "\ndaemon off;" >> /etc/nginx/nginx.conf
RUN sed -i -e "s/;\?daemonize\s*=\s*yes/daemonize = no/g" /etc/php/7.2/fpm/php-fpm.conf

RUN apt-get install -y dialog apt-utils
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:inkscape.dev/stable
RUN apt update
RUN apt-get -y install inkscape
RUN apt-get -y install wget
RUN apt-get -y install unzip

RUN apt-get -y install imagemagick && \
  apt-get -y install nodejs && \
  apt-get -y install npm

WORKDIR /var/www/html
COPY . /var/www/html

# Nginx config
RUN rm -rf /etc/nginx/sites-available && rm -rf /etc/nginx/sites-enabled/default
COPY default /etc/nginx/sites-enabled/default

RUN npm install -g bower && \
  npm install -g grunt-cli && \
  npm install grunt --save-dev && \
  npm install grunt-contrib-clean grunt-contrib-uglify grunt-text-replace grunt-contrib-less grunt-svgmin grunt-contrib-copy grunt-contrib-watch

# Make our shell script executable
RUN chmod +x env.sh

RUN ./env.sh
RUN echo '{ "allow_root": true }' > /root/.bowerrc
RUN chmod a+x build.sh
RUN ./build.sh fetch jmol

# CMD ["/bin/bash", "-c", "/var/www/html/env.sh && /etc/init.d/php7.2-fpm start && nginx -t"]
CMD ["/bin/bash", "-c", "/var/www/html/env.sh && /etc/init.d/php7.2-fpm start && /etc/init.d/nginx start"]

# Expose ports.
EXPOSE 80