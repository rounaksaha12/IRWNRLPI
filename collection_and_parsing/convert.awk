#!/usr/bin/gawk -f

#
#   run as:
#   ./extract.awk combined_prot_table.tsv trimmed_combined_table.tsv ENCODE.tsv 
#

function fasta(_name,_idx,_sequence,_filename, _i){
    print ">"_idx"\t"_name > _filename;
    _i=1;
    while(_i<=length(_sequence)){
        print substr(_sequence,_i,50) > _filename;
        _i+=50;
    }
}

BEGIN{
    FS="\t";
}

FILENAME=="combined_prot_table.tsv" && FNR>1{
    split($2,names,",");
    for(name in names){
        id[names[name]]=$3;
        prot_seq[$3]=$4;
    }
}

FILENAME=="trimmed_combined_table.tsv" && FNR>1{
    split($4,nonhsat,".");
    split($7,nonhsag,".");

    if(nonhsat[1] in known_rna){}
    else known_rna[nonhsat[1]]=$8;

    if(nonhsag[1] in known_rna){}
    else known_rna[nonhsag[1]]=$8;
}

FILENAME=="sampled_rna.txt"{
    fasta($2,$1,known_rna[$2],"rna_list_with_sequences.txt");
    rna_identifier[$2]=$1;
}

FILENAME=="sampled_prot.txt"{
    fasta($2,$1,prot_seq[$2],"prot_list_with_sequences.txt");
    prot_identifier[$2]=$1;
}

FILENAME=="sampled_inter.tsv" && FNR>1{
    print "hello"
    print rna_identifier[$3]"\t"$3"\t"prot_identifier[$5]"\t"$5 > "list_of_interactions.txt";
}