for a in \
CY006079 \
CY005970 \
; do
    /net/intdev/oblast01/infernal/notebook/23_0925_vadr_1p6_release/test-install/infernal/binaries/cmbuild -n $a --verbose --noss --hand --noh3pri $a.vadr.cm $a.2.stk
done


