#!/bin/bash
../../examples/generate-pattern.x collections/KERNEL256/text.TEXT 200 9 1 > KERNELpat.txt
../../examples/generate-pattern.x collections/CC256/text.TEXT 200 9 1 > CCpat.txt
cp KERNELpat.txt collections/KERNEL16/patterns/words.KERNEL16.9.words.txt
cp KERNELpat.txt collections/KERNEL64/patterns/words.KERNEL64.9.words.txt
cp KERNELpat.txt collections/KERNEL256/patterns/words.KERNEL256.9.words.txt
cp KERNELpat.txt collections/KERNEL1024/patterns/words.KERNEL1024.9.words.txt
cp KERNELpat.txt collections/KERNEL4096/patterns/words.KERNEL4096.9.words.txt
cp KERNELpat.txt collections/KERNEL16384/patterns/words.KERNEL16384.9.words.txt
cp CCpat.txt collections/CC16/patterns/words.CC16.9.words.txt
cp CCpat.txt collections/CC64/patterns/words.CC64.9.words.txt
cp CCpat.txt collections/CC256/patterns/words.CC256.9.words.txt
cp CCpat.txt collections/CC1024/patterns/words.CC1024.9.words.txt
cp CCpat.txt collections/CC4096/patterns/words.CC4096.9.words.txt
cp CCpat.txt collections/CC16384/patterns/words.CC16384.9.words.txt
