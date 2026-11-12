#!/bin/bash
# Consumes the producer's output; fails loudly when it is missing.
set -e
wc -l < out.txt > result.txt
