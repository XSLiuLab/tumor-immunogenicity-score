#!/bin/sh

cp -f report/comp_profile.png report/report_main.html docs/
cp -rf report/report_main_files/ docs/report_main_files
mv docs/report_main.html docs/index.html

