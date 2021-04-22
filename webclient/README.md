# MolView

This repository holds the code of MolView.org.
It is not actively maintained due to lacking financial support.

## Running webclient dev

0. use CMD ["/bin/bash", "-c", " /etc/init.d/php7.2-fpm start && /etc/init.d/nginx start"] in dockerfile.
1. docker build -t molsearch .
2. cd webclient
3. docker run -td -p 8001:80 -e API_URL=http://192.168.1.85:35001 -v "$PWD":/var/www/html molsearch
4. npx ./build.sh (build app)

## License

```
MolView (http://molview.org)
Copyright (c) 2014, 2015 Herman Bergwerf

MolView is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MolView is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with MolView.  If not, see <http://www.gnu.org/licenses/>.
```

See LICENSE.md for more details about third party code.

## Build

Instructions:

1. Copy all files to a server with support for: PHP5, php5-curl, .htaccess and mod_rewrite
2. Configure $directory in `php/utility.php` to point to the MolView root directory
3. Configure ErrorDocument in `.htaccess` to point to `page.php`
4. Make sure the Inkscape and the ImageMagick CLI are installed
5. Install npm and install local npm modules
6. Run the build script (use `./build.sh fetch jmol` for a clean build)
   - **Only run bower and grunt, and render images:** `./build.sh`
   - **Also fetch external PHP sources:** `./build.sh fetch`
   - **Also fetch Jmol from sourceforge:** `./build.sh fetch jmol`
   - **Also fetch Jmol from stolaf.edu:** `./build.sh fetch jmol nightly`
   - **Only render images:** `./build.sh render`
