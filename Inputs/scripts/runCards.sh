#!/bin/bash

YEAR=$1

./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 0 mupi
./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 0 murho
./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 0 mua1

./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 1 mupi
./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 1 murho
./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 1 mua1

./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 2 mupi
./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 2 murho
./HTC_submit_cards.sh card_config_newDNN_$YEAR.conf 2 mua1

