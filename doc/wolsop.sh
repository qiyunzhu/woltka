#!/bin/bash

#######################################
#
# Standard workflow for WoL + Woltka
#
#######################################

# Author:  Qiyun Zhu
# License: BSD-3-Clause
# Version: 0.0.1-dev
# Email: qiyunzhu@gmail.com

# Last updated: 2021-04-17

# Usage: Customize the "Parameters" section, then run this script.


##############
# Parameters #
##############

# WoL data release directory
#   available for download from the WoL Globus endpoint
#   see WoL doc for instruction
db=

# input alignment file or directory
#   see Woltka doc for instruction
input=

# input file extension (optional)
ext=

# output format (biom or tsv)
fmt=biom

# taxonomy tree format (taxdump or lineage)
taxtree=taxdump

# taxonomic ranks
taxranks=phylum,class,order,family,genus,species

# GO classification (yes/no)
go=yes

# MetaCyc classification (yes/no)
metacyc=yes

# KEGG classification (specify KEGG directory path)
#   this directory is not part of WoL data release, see Woltka doc for
#   instruction
kegg=


#######
# SOP #
#######

[ "$fmt" == biom ] && altfmt= || altfmt="--to-tsv"
[ -z "$ext" ] && filext= || filext="--filext $ext"

# taxonomic classification (NCBI)

if [ "$taxtree" == taxdump ]; then
  woltka classify \
    --input  $input \
    --map    $db/taxonomy/taxid.map \
    --nodes  $db/taxonomy/nodes.dmp \
    --names  $db/taxonomy/names.dmp \
    --rank   none,free,$taxranks \
    $filext \
    $altfmt \
    --output .

else
  woltka classify \
    --input   $input \
    --lineage $db/taxonomy/lineage.txt \
    --rank    none,free,$taxranks \
    $filext \
    $altfmt \
    --output .
fi

mv none.$fmt ogu.$fmt

# functional classification (UniRef)

woltka classify \
  --input  $input \
  --coords $db/proteins/coords.txt.xz \
  --map    $db/function/uniref/uniref.map.xz \
  --names  $db/function/uniref/uniref.name.xz \
  --map-as-rank \
  --rank   none,uniref \
  $filext \
  $altfmt \
  --output .

mv none.$fmt orf.$fmt


######
# GO #
######

if [ "$go" == yes ]; then

  go=$db/function/go

  mkdir -p go
  cd go

  # by domain
  for domain in all component function process; do
    woltka tools collapse -m $go/$domain.map.xz -n $go/name.txt \
      -i ../uniref.$fmt -o $domain.$fmt

    # GO slim (generic)
    woltka tools collapse -m $go/generic/$domain.map -n $go/name.txt \
      -i $domain.$fmt -o $domain.generic.$fmt
  done
  cd ..
fi


###########
# MetaCyc #
###########

if [ "$metacyc" == yes ]; then

  mc=$db/function/metacyc

  mkdir -p metacyc
  cd metacyc

  # ORF to protein
  woltka tools collapse -m $mc/protein.map.xz -n $mc/protein_name.txt \
    -i ../orf.$fmt -o protein.$fmt

  # protein to enzrxn (enzymatic reaction)
  woltka tools collapse -m $mc/protein-to-enzrxn.txt -n $mc/enzrxn_name.txt \
    -i protein.$fmt -o enzrxn.$fmt

  # enzrxn to reaction
  woltka tools collapse -m $mc/enzrxn-to-reaction.txt -n $mc/reaction_name.txt \
    -i enzrxn.$fmt -o reaction.$fmt

  # reaction to pathway
  woltka tools collapse -m $mc/reaction-to-pathway.txt -n $mc/pathway_name.txt \
    -i reaction.$fmt -o pathway.$fmt

  # pathway to super pathway
  woltka tools collapse -m $mc/pathway-to-super_pathway.txt -n $mc/pathway_name.txt \
    -i pathway.$fmt -o super_pathway.$fmt

  # super pathway (or pathway) to pathway type
  woltka tools collapse -m $mc/pathway_type.txt -n $mc/all_class_name.txt \
  -i super_pathway.$fmt -o pathway_type.$fmt

  # protein to gene
  woltka tools collapse -m $mc/protein-to-gene.txt -n $mc/gene_name.txt \
    -i protein.$fmt -o gene.$fmt

  # protein to go
  woltka tools collapse -m $mc/protein-to-go.txt \
    -i protein.$fmt -o go.$fmt

  # enzrxn to regulation to regulator
  woltka tools collapse -m $mc/enzrxn-to-regulation.txt \
    -i enzrxn.$fmt -o regulation.$fmt
  woltka tools collapse -m $mc/regulation-to-regulator.txt -n $mc/compound_name.txt \
    -i regulation.$fmt -o regulator.$fmt

  # reaction to compound (left and right)
  for side in left right; do
    woltka tools collapse -m $mc/reaction-to-${side}_compound.txt -n $mc/compound_name.txt \
      -i reaction.$fmt -o ${side}_compound.$fmt

    # compound type
    woltka tools collapse -m $mc/compound_type.txt -n $mc/all_class_name.txt \
      -i ${side}_compound.$fmt -o ${side}_compound_type.$fmt
  done

  # compound and type (both sides)
  woltka tools merge -i left_compound.$fmt -i right_compound.$fmt -o compound.$fmt
  woltka tools merge -i left_compound_type.$fmt -i right_compound_type.$fmt -o compound_type.$fmt

  # reaction to EC
  woltka tools collapse -m $mc/reaction-to-ec.txt \
    -i reaction.$fmt -o ec.$fmt

  # pathway coverage (by reaction)
  woltka tools coverage -m $mc/pathway-to-reaction_list.txt \
    -i reaction.$fmt -o pathway_coverage.$fmt

  cd ..

fi


########
# KEGG #
########

if [ -d "$kegg" ]; then

  ke=$kegg

  mkdir -p kegg
  cd kegg

  # entrance
  # UniRef to KO
  woltka tools collapse -m $db/function/kegg/ko.map.xz -n $db/function/kegg/ko.name \
    -i ../uniref.$fmt -o ko.$fmt

  # main cascade
  # KO to reaction
  woltka tools collapse -m $ke/ko-to-reaction.txt -n $ke/reaction_name.txt \
    -i ko.$fmt -o reaction.$fmt

  # reaction to module
  woltka tools collapse -m $ke/reaction-to-module.txt -n $ke/module_name.txt \
    -i reaction.$fmt -o module.$fmt

  # module to pathway
  woltka tools collapse -m $ke/module-to-pathway.txt -n $ke/pathway_name.txt \
    -i module.$fmt -o pathway.$fmt

  # classes
  # reaction to rclass
  woltka tools collapse -m $ke/reaction-to-rclass.txt -n $ke/rclass_name.txt \
    -i reaction.$fmt -o rclass.$fmt

  # module class
  woltka tools collapse -m $ke/module-to-class.txt \
    -i module.$fmt -o module_class.$fmt

  # pathway class
  woltka tools collapse -m $ke/pathway-to-class.txt \
    -i pathway.$fmt -o pathway_class.$fmt

  # compounds
  for side in left right; do
    woltka tools collapse -m $ke/reaction-to-${side}_compound.txt -n $ke/compound_name.txt \
      -i reaction.$fmt -o ${side}_compound.$fmt
  done
  woltka tools merge -i left_compound.$fmt -i right_compound.$fmt -o compound.$fmt

  # extended
  # KO to EC
  woltka tools collapse -m $ke/ko-to-ec.txt \
    -i ko.$fmt -o ec.$fmt

  # KO to GO
  woltka tools collapse -m $ke/ko-to-go.txt \
    -i ko.$fmt -o go.$fmt

  # KO to COG
  woltka tools collapse -m $ke/ko-to-cog.txt \
    -i ko.$fmt -o cog.$fmt

  # KO to disease
  woltka tools collapse -m $ke/ko-to-disease.txt -n $ke/disease_name.txt \
    -i ko.$fmt -o disease.$fmt

  # coverage
  # module coverage by reaction
  woltka tools coverage -m $ke/module-to-reaction.txt \
    -i reaction.$fmt -o module_coverage.$fmt

  # pathway coverage by module
  woltka tools coverage -m $ke/pathway-to-module.txt \
    -i module.$fmt -o pathway_coverage.$fmt

  cd ..

fi
