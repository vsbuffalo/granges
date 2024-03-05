#!/bin/bash

START_SEED=1
END_SEED=10000

for SEED in $(seq $START_SEED $END_SEED); do
  echo "Running tests with TEST_SEED=$SEED"
  TEST_SEED=$SEED cargo test --test bedtools_validation
  if [ $? -ne 0 ]; then
    echo "Tests failed with TEST_SEED=$SEED"
    exit 1
  else
    echo "Tests passed with TEST_SEED=$SEED"
  fi
done

