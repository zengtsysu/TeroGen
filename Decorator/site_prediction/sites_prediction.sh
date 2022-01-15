#!/bin/sh

echo " Job begin time: " `date`

### build vocab for training and validation set
#onmt_build_vocab -config vocab.yaml -n_sample -1

### train the model
#onmt_train -config config.yaml

### make prediction
onmt_translate -model data/checkpoints/sites_can.pt data/checkpoints/sites_mix.pt\
               -src data/decoration/test_src.txt -output data/decoration/result/pred_test_ensemble.txt \
               -batch_size 64 -replace_unk -max_length 300 -beam_size 10 -n_best 10 -gpu 0 \

echo " Job finish time: " `date`
