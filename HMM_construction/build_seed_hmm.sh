#!/bin/bash

muscle5 \
    -align cloned_kfp_beta_fingers.faa \
    -output cloned_kfp_beta_fingers.aln

hmmbuild \
    beta_finger.seed.hmm \
    cloned_kfp_beta_fingers.aln 