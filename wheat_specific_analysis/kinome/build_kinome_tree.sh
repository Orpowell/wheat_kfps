#!/bin/bash

muscle5 \
	-super5 wheat_chinese_spring_kinome.faa \
	-output wheat_chinese_spring_kinome.aln

fasttree wheat_chinese_spring_kinome.aln > wheat_chinese_spring_kinome.nwk
