#!/bin/bash
# Get the current date and time
current_date_time=$(date +"%Y-%m-%d %H:%M:%S")
git add .
git commit -m "Modification at $current_date_time."
git push -u origin main --force
