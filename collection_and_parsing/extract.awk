#!/usr/bin/gawk -f

#
#   run as:
#   ./extract.awk combined_prot_table.tsv trimmed_combined_table.tsv HT_without_pred.tsv 
#

function fasta(_name,_idx,_sequence,_filename, _i){
    print ">"_idx" "_name > _filename;
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

FILENAME=="HT_without_pred.tsv" && $7=="Homo sapiens" && $4=="lncRNA" && $6=="protein"{
    if($3=="-" || $5=="-") next;

    if(!($3 in known_rna)) next;

    split($5,prot_arr,";");
    
    for(protein in prot_arr){
        if(!(prot_arr[protein] in id)) continue;
        output_string=$3"\t"id[prot_arr[protein]];
        if(output_string in unique_inter_pair){
            continue;
        }
        unique_inter_pair[output_string]=1;
        if($3 in rna_list){
            rna_list[$3]++;
        }
        else{
            rna_list[$3]=1;
        }

        if(id[prot_arr[protein]] in prot_list){
            prot_list[id[prot_arr[protein]]]++;
        }
        else{
            prot_list[id[prot_arr[protein]]]=1;
        }

    }

}

END{
    print length(known_rna)
    print length(id)
    print length(prot_list)
    rna_cnt=0;
    prot_cnt=0;

    for(pair in unique_inter_pair){
        split(pair,out,"\t");
        if(rna_list[out[1]]>1 && prot_list[out[2]]>1) continue;
        if(rna_list[out[1]]==1 && prot_list[out[2]]==1) continue;
        #print "halo "out[1]": "rna_list[out[1]]", "out[2]": "prot_list[out[2]];
        if(rna_list[out[1]]>1 && prot_list[out[2]]<=1) rna_list[out[1]]--;
        else if(rna_list[out[1]]<=1 && prot_list[out[2]]>1) prot_list[out[2]]--;
    }

    for(rna in rna_list){
        if(rna_list[rna]>1){
            #print rna_cnt"\t"rna"\t"known_rna[rna] > "rounak_rna.txt";
            fasta(rna,rna_cnt,known_rna[rna],"CTAON_rna.txt");
            rna_mapping[rna]=rna_cnt;
            rna_cnt+=1;
        }
    }

    for(prot in prot_list){
        if(prot_list[prot]>1){
            #print prot_cnt"\t"prot"\t"prot_seq[prot] > "rounak_prot.txt";
            fasta(prot,prot_cnt,prot_seq[prot],"CTAON_prot.txt");
            prot_mapping[prot]=prot_cnt;
            prot_cnt+=1;
        }
    }

    for(pair in unique_inter_pair){
        split(pair,out,"\t");
        if(rna_list[out[1]]<=1 || prot_list[out[2]]<=1) continue;
        #print rna_mapping[out[1]]"\t"out[1]"\t"prot_mapping[out[2]]"\t"out[2] > "rounak_interactions.txt"
        print rna_mapping[out[1]]"\t"out[1]"\t"prot_mapping[out[2]]"\t"out[2] > "CTAON_interactions.txt"
    }
}
