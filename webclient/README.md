Build
---
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
