#for a in GM968017; do
#for a in GM968018; do
#for a in GM968019; do
#for a in M35997; do 
for a in GM968015; do
    esl-sfetch --index $a.vadr.protein.fa 
    for c in 1..728:+,957..957:+; do 
        esl-sfetch $a.vadr.protein.fa $a.1/$c | $VADRHMMERDIR/hmmbuild -n $a/$c --informat afa tmp.1.hmm -
    done
    for c in 706..1125:+; do 
        esl-sfetch $a.vadr.protein.fa $a.1/$c | $VADRHMMERDIR/hmmbuild -n $a/$c --informat afa tmp.2.hmm -
    done
done
