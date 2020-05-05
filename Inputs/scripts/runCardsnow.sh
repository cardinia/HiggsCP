p#!/bin/bash

YEAR=$1

./NohupCardsParser.sh $YEAR 0 mupi 
./NohupCardsParser.sh $YEAR 0 murho 
./NohupCardsParser.sh $YEAR 0 mua1 
./NohupCardsParser.sh $YEAR 0 mu0a1 

./NohupCardsParser.sh $YEAR 0 
./NohupCardsParser.sh $YEAR 1 
./NohupCardsParser.sh $YEAR 2 


./NohupCardsParser.sh $YEAR 1 mupi 
./NohupCardsParser.sh $YEAR 1 murho 
./NohupCardsParser.sh $YEAR 1 mua1 
./NohupCardsParser.sh $YEAR 1 mu0a1 

./NohupCardsParser.sh $YEAR 2 mupi 
./NohupCardsParser.sh $YEAR 2 murho 
./NohupCardsParser.sh $YEAR 2 mua1 
./NohupCardsParser.sh $YEAR 2 mu0a1 
