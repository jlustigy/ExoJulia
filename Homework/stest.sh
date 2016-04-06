#!/bin/bash
if [ $# -ne 1 ]
then
    echo "Please provide the homework folder to test (e.g. hw1)"
    exit
fi

echo ""
echo "Calculating Speeds for: $1"

for HW in `find . -type d`
do
   if [ "$HW" == "./$1" ]
   then
      for GROUP in `find $HW -type d`
      do
         if [ "$GROUP" != "$HW" ]
         then
            echo ""
            echo "-----------------------"
            echo "Running: $GROUP"
            for FILE in `find $GROUP -name "*.jl"`
            do
               sed -n 's/^#@stest[[:space:]]*//p' "$FILE" | while read -r line ; do
                  echo ""
                  echo "  Runtime for: $line"
                  echo "  $(julia -e"include(\"$FILE\"); @time $line")"
               done
            done
         fi
      done
   fi
done
echo ""
