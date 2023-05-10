$file=$ARGV[0];    # input a fasta file of genome/metagenome
################ 
## 1. setup “prodigal.linux” , "bedtools" and “pilercr” in the default environment
## 2. example:  perl  0.Cas-Finder.pl  $sample.fasta 
################ 


$number="$file\_metacontig\_";
################
$filein="$file";
$fileout="$file.3kb.fa";
open(FASTA,"$filein");
open(FA,">$fileout");
$s=0;
while($line=<FASTA>)
{  chomp $line;
  if($line=~/^>/)
   { if($s>0)
     { if(length($seq)>3000){print FA "$pre\n$seq\n";} }
       $pre=$line;  $seq="";
   }
    else{$seq="$seq$line";}
	$s++;
}
if(length($seq)>3000)
{print FA "$pre\n$seq\n";}
close FASTA; close FA;
######################
open(FF,"$file.3kb.fa");
open(F1,">$file.crisper.scaffold.fa");  close F1;
open(F2,">$file.pep.fasta");       close F2;
open(F3,">$file.pep.cas.fasta");   close F3;
open(F4,">$file.crispr.loc");      close F4;
open(F5,">$file.crispr.spacer");   close F5;
open(AAA,">temp");
$n=0;  $ss=0;
while($myline=<FF>)
{  chomp $myline;  $n++;  
   if($myline=~/^>/)
   { 
     if(length($seq)>3000  && $ss>0){print AAA "$head\n$seq\n";   }   
	 $head=$myline;  $seq="";   $ss=1;
	 ############################################################################################################################################
	 if($n>10000)
	 {
	  close AAA;
     ############################################## 1.cas & protein annotation
      `prodigal.linux   -a  temp.pep   -i  temp  -p single  -f gff  -o  temp.gff`;
      `pilercr  -in  temp  -out   temp.spacer`;  
     ########################  2.get cas locs
      open(A,"temp.spacer");
      open(B,">temp.spacer.loc");
      while($line=<A>)
      { chomp $line;
       if($line=~/^SUMMARY BY POSITION/){$s=1;}
       if($line=~/^Help on reading this report/){$s=0;}
       if($s>0)
        {  
          @bb=split/\s+/,$line;  $len=@bb;
          if($line=~/^>/){ @ar=split/\s+/,$line; $scaffold=$ar[0];  $scaffold=~s/>//;}
          if($len>7 && $bb[1]>0 )
             {  
              if($start1<0){$start1=1;}
              if($pre eq "=")
                { $start=$bb[$len-6];  
                  $end=$bb[$len-6]+$bb[$len-5]+1;
                  $start1=$start-10000;
                  $end1=$end+10000;
                  if($start1<1){$start1=1;}
                  print B "$number\_",$scaffold,"\t",$start1,"\t",$end1,"\t",$start,"\t",$end,"\t",$bb[$len-1],"\t",$bb[$len-4],"\t",$bb[$len-3],"\t",$bb[$len-2],"\t",$bb[1],"\n";
                }
            else{ $start=$bb[$len-7];
                  $end=$bb[$len-7]+$bb[$len-6]+1;
                  $start1=$start-10000;
                  $end1=$end+10000;
                  if($start1<1){$start1=1;}
                 print B "$number\_",$scaffold,"\t",$start1,"\t",$end1,"\t",$start,"\t",$end,"\t",$bb[$len-1],"\t",$bb[$len-5],"\t",$bb[$len-4],"\t",$bb[$len-3],"\t",$bb[1],"\n";
                }
             }
          @pp=split//,$line;  $pre=$pp[0];  
        }
      }
      close A;   close B; 
	 
      #######################  3.get protein locs
      open(A,"temp.gff");
      open(B,">temp.gff.loc");
      while($line=<A>)
      { chomp $line;  
        @bb=split/\t/,$line;
        if($bb[2] eq "CDS")
        {  $scaffold="$number\_$bb[0]";
           $ll=$bb[8];  $ll=~s/;/\t/g;  $ll=~s/\_/\t/g;
           @b=split/\t/,$ll; 
           $gene="$scaffold\_$b[1]";
           print B "$scaffold\t$bb[3]\t$bb[4]\t$bb[6]\t$gene\n";
           $loc{$gene}="$scaffold\t$bb[3]\t$bb[4]\t$bb[6]";
        }
      }
      close A; close B; 
      ######################################  4.cas locs vs protein locs
      `bedtools intersect -wo -a temp.gff.loc    -b  temp.spacer.loc  >  temp.bed`;
      open(A,"temp.bed");
      while($line=<A>)
      {  chomp $line;
         @bb=split/\s+/,$line;  
         $line=~s/\_metacontig\_\_/\t/g;  
         @ar=split/\s+/,$line;
         $good{$bb[4]}=1;    
         $cpscaffold{$ar[1]}=1;   #print $ar[1],"\n";
      }
      close A;
      ######################################  5.get cas locs protein
      open(A,"temp.pep");
      open(B,">temp.pep.fasta");
      open(C,">temp.pep.filtered.fasta");
      while($line=<A>)
      { chomp $line;
        $line=~s/>/>$number\_/;
        if($line=~/^>/)
         { $ll=$line; $ll=~s/>//;  
           @bb=split/\s+/,$ll; #print $bb[0],"\n";
           print B ">",$bb[0],"\t",$loc{$bb[0]},"\n";
           if(exists($good{$bb[0]})){$s=1; print C ">",$bb[0],"\t",$loc{$bb[0]},"\n";}
           else{$s=0;}
         }
         else{print B $line,"\n";
              if($s>0){print C $line,"\n";}
             } 
      }
      close A; close B; close C;

	  ####################################### 6.get cas locas scaffold
                  open(A,"temp");
                  open(B,">temp-crisper-scaffold.fa");
                  while($line=<A>)
                   { $ll=$line; chomp $ll; $ll=~s/>//;  @bb=split/\s+/,$ll;
                     $line2=<A>;
                    if(exists($cpscaffold{$bb[0]})){print B $line,$line2;}
                   }
                  close A; close B;

	  ####################################### 7.save result files
	  `cat  temp-crisper-scaffold.fa  >> $file.crisper.scaffold.fa`;
	  `cat  temp.pep.fasta   >>  $file.pep.fasta`;
              `cat  temp.pep.filtered.fasta  >>  $file.pep.cas.fasta`;
	  `cat  temp.spacer.loc  >>  $file.crispr.loc`;
              `cat  temp.spacer      >>  $file.crispr.spacer`;
	  `rm temp*`;
              `rm  $file `;  
	  #############################################################################################################################
	  $n=0;  
	  open(AAA,">temp");
	 }
   }
   else{$seq="$seq$myline";}
}
close FF;  

##########
if(length($seq)>3000  && $ss>0){print AAA "$head\n$seq\n"; }
close AAA;
############################################## 1.cas & protein annotation
`prodigal.linux   -a  temp.pep   -i  temp  -p single  -f gff  -o  temp.gff`;
`pilercr  -in  temp  -out   temp.spacer`;
########################  2.get cas locs
open(A,"temp.spacer");
open(B,">temp.spacer.loc");
while($line=<A>)
{ chomp $line;
  if($line=~/^SUMMARY BY POSITION/){$s=1;}
  if($line=~/^Help on reading this report/){$s=0;}
  if($s>0)
    {  
     @bb=split/\s+/,$line;  $len=@bb;
     if($line=~/^>/){ @ar=split/\s+/,$line; $scaffold=$ar[0];  $scaffold=~s/>//;}
     if($len>7 && $bb[1]>0 )
      {  
        if($start1<0){$start1=1;}
        if($pre eq "=")
         { $start=$bb[$len-6];  
           $end=$bb[$len-6]+$bb[$len-5]+1;
           $start1=$start-10000;
           $end1=$end+10000;
           if($start1<1){$start1=1;}
           print B "$number\_",$scaffold,"\t",$start1,"\t",$end1,"\t",$start,"\t",$end,"\t",$bb[$len-1],"\t",$bb[$len-4],"\t",$bb[$len-3],"\t",$bb[$len-2],"\t",$bb[1],"\n";
         }
        else{ $start=$bb[$len-7];
              $end=$bb[$len-7]+$bb[$len-6]+1;
              $start1=$start-10000;
              $end1=$end+10000;
            if($start1<1){$start1=1;}
              print B "$number\_",$scaffold,"\t",$start1,"\t",$end1,"\t",$start,"\t",$end,"\t",$bb[$len-1],"\t",$bb[$len-5],"\t",$bb[$len-4],"\t",$bb[$len-3],"\t",$bb[1],"\n";
            }
        }
       @pp=split//,$line;  $pre=$pp[0];  
     }
}
close A;   close B; 
	 
#######################  3.get protein locs
open(A,"temp.gff");
open(B,">temp.gff.loc");
while($line=<A>)
{ chomp $line;  
  @bb=split/\t/,$line;
  if($bb[2] eq "CDS")
   {  $scaffold="$number\_$bb[0]";
      $ll=$bb[8];  $ll=~s/;/\t/g;  $ll=~s/\_/\t/g;
      @b=split/\t/,$ll; 
      $gene="$scaffold\_$b[1]";
      print B "$scaffold\t$bb[3]\t$bb[4]\t$bb[6]\t$gene\n";
      $loc{$gene}="$scaffold\t$bb[3]\t$bb[4]\t$bb[6]";
   }
}
close A; close B; 
######################################  4.cas locs vs protein loc
`bedtools intersect -wo -a temp.gff.loc    -b  temp.spacer.loc  >  temp.bed`;
open(A,"temp.bed");
while($line=<A>)
{  chomp $line;
   @bb=split/\s+/,$line;  
   $line=~s/\_metacontig\_\_/\t/g;  
   @ar=split/\s+/,$line;
   $good{$bb[4]}=1;    
   $cpscaffold{$ar[1]}=1; #print $ar[1],"\n";
}
close A;
######################################  5.get cas locs protein
open(A,"temp.pep");
open(B,">temp.pep.fasta");
open(C,">temp.pep.filtered.fasta");
while($line=<A>)
{ chomp $line;
  $line=~s/>/>$number\_/;
  if($line=~/^>/)
    { $ll=$line; $ll=~s/>//;  
      @bb=split/\s+/,$ll; #print $bb[0],"\n";
      print B ">",$bb[0],"\t",$loc{$bb[0]},"\n";
      if(exists($good{$bb[0]})){$s=1; print C ">",$bb[0],"\t",$loc{$bb[0]},"\n";}
      else{$s=0;}
     }
   else{ print B $line,"\n";
         if($s>0){print C $line,"\n";}
       } 
 }
close A; close B; close C;
####################################### 6.get cas locas scaffold
open(A,"temp");
open(B,">temp-crisper-scaffold.fa");
while($line=<A>)
{ $ll=$line; chomp $ll; $ll=~s/>//;  @bb=split/\s+/,$ll;
  $line2=<A>;
  if(exists($cpscaffold{$bb[0]})){print B $line,$line2;}
}
close A; close B;
####################################### 7.save result files
`cat  temp-crisper-scaffold.fa  >> $file.crisper.scaffold.fa`;
`cat  temp.pep.fasta   >>  $file.pep.fasta`;
`cat  temp.pep.filtered.fasta  >>  $file.pep.cas.fasta`;
`cat  temp.spacer.loc  >>  $file.crispr.loc`;
`cat  temp.spacer      >>  $file.crispr.spacer`;   
`rm temp*`;
#`rm  $file `;  
#################################################################################################################################
