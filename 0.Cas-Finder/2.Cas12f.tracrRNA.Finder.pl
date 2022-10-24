#Example: perl 2.Cas12f.tracrRNA.Finder.pl  Contig_name   Cas12f_start  Cas12f_end  Direction CRISPR_array_start   CRISPR_array_end   Direct_repeat    output_file  

$reference=$ARGV[0];    # the fasta file of  genome
$contig=$ARGV[1];
$gene_start=$ARGV[2];
$gene_end=$ARGV[3];
$direction=$ARGV[4];
$array_start=$ARGV[5];
$array_end=$ARGV[6];
$Direct_repeat=$ARGV[7];
$output=$ARGV[8];
if($output eq ""){$outfile="$contig.tracrRNA.txt";}
else{$outfile="$output.tracrRNA.txt";}
######################
open(FASTA,"$reference");
open(FA,">$reference.fa");
$s=0;
while($line=<FASTA>)
{ chomp $line;
  if($line=~/^>/)
   { if($s>0){print FA "\n";}
      print FA $line,"\n";
     $s=1;
   }
    else{print FA $line;}
}    print FA "\n";
close  FASTA; close FA;
#####
open(AA,"$reference.fa");
while($line=<AA>)
{  $line2=<AA>;
   chomp $line;  chomp $line2;  
   $ll=$line;  $ll=~s/>//;
   @bb=split/\s+/,$ll;
   $h{$bb[0]}=$line2;
}
close AA;
`rm  $reference.fa`;
######################
if($array_start>$gene_start){  $p1=$gene_start;  $p2=$gene_end;  $p3=$array_start;  $p4=$array_end; }
else{ $p3=$gene_start;  $p4=$gene_end;  $p1=$array_start;  $p2=$array_end; }
$lenn=$p3-$p2+1;     
$line=substr($h{$contig},$p2,$lenn); 
if($direction eq "-"){   $line=~tr/ATCGatcg-/TAGCTAGC-/;    $line=reverse($line);  }
##################################################
$line2=$Direct_repeat;  chomp $line2;   
$number=22;  
$ii=length($line2)-$number;   
$grna=substr($line2,$ii,$number);
#print $grna,"\n";

$sca=200; #### extend length


$alen=length($line);
$ii=0;
open(BB,">$outfile");
while($ii<$alen)
{ 
      $seq=substr($line,$ii,$number);
      $seq2=$seq;   $seq2=~tr/ATCGatcg/TAGCTAGC/; $seq2=reverse($seq2);
      #$seqq=substr($line,$ii,$sca);
      if($ii-$sca+$number>=0){$start=$ii-$sca+$number;  $seqq=substr($line,$start,$sca);}
       else{$start=1;  $seqq=substr($line,$start,$ii-$sca+$number+$sca-1);   }
      $seqq2=$seqq;     # $seqq2=~tr/ATCGatcg/TAGCTAGC/;   $seqq2=reverse($seqq2);
      ###
       #$seqqf=substr($line,$ii-$sca+$number,$sca);  
       $seqqf=substr($line,$ii,$sca);
       $seqqf=~tr/ATCGatcg/TAGCTAGC/; $seqqf=reverse($seqqf);
      #$ff1=getsite($seq,$ii,'F',$grna,$seqqf,$contig);   #### Forword
      $ff2=getsite($seq2,$ii,'R',$grna,$seqq2,$contig);   #### Reverse
       #print BB $ff1,$ff2;
      @ar=split/\n/,$ff2;   @ll=split/\t/,$ar[0];   @ll1=split//,$ar[3];  @ll2=split//,$ar[2];  @ll3=split//,$ar[1];
      $t=0;  $i=0;  $bad=0; $len=@ll1;  $dr="";  $match="";  $good=0;  $tr=""; 
      while($i<=$len)
      {    $temp="$ll1[$i]$ll1[$i+1]$ll1[$i+2]";
           if($temp eq "***"){$t=1;}
           if($t>0){$dr="$dr$ll2[$i]";   $match="$match$ll1[$i]";  $tr="$tr$ll3[$i]";
                              if($ll2[$i] eq $ll3[$i]){$good++; }
                   }
                        elsif($ll3[$i] ne "-"){$bad++;}
            $i++;
      }         
                   $trancrRNA=substr($ar[4],0,length($ar[4])-2);
                   $match2=$match; $match2=~s/\ \*\ //g; $match2=~s/\ \*\*\ //g;  $match2=~s/\ \*\*\ \ //g;     $match2=~s/\ \ \*\*\ //g;    $match2=~s/\ //g;  $rate=(length($match2)+1)/(length($match)+1);
                     $tr2=substr($tr,0,10);   
                    if($t>0  &&  $good>=9  &&  $rate>0.65  &&  length($trancrRNA)>130    &&   not  exists($v{$tr2}))
                               { $site=$p2+$ii;     
                  print  BB  "Location   $contig   $site\t$ii\n";
                  print BB  "Anti-repeat     $tr\nDirect repeat:  $dr\nMatch:          $match\nGood match:  $good\nTrancrRNA:   $trancrRNA\n";  
                  print BB "Scaffold:    $trancrRNA  GAAA  $dr\n\n\n";
                          $v{$tr2}=1;   
                                   # if($ff2){$ii=$ii+30;}
               }     
      $ii++;
}
close BB;


sub  getsite
{
 ($find,$loc,$dir,$gRNA,$longseq,$con)=@_;
 @ar=split//,$find;
 @br=split//,$gRNA;
 $len=length($gRNA);
 $i=0;  $maxp=0;  $good="";
 while($i<$len)
 {
   $a1=substr($gRNA,0,$i);
   $b1=substr($gRNA,$i,$len-$i+1);
   $a2=substr($find,0,$i);
   $b2=substr($find,$i,$len-$i+1);
   $t1="-";  $t2="--";
############  gRNA with one insertion
   $sum1="$a1$t1$b1";
   @temp1=split//,$sum1;
   $j=0;  $ss=""; $p=0;
   while($j<($len+1))
    {    $tt="$ar[$j]$temp1[$j]";
         if($ar[$j] eq $temp1[$j]  &&  $ar[$j] ne ""){$ss="$ss*";  $p++;}
          elsif($tt eq "AG" || $tt eq "CT"){$ss="$ss*";  $p++;}
          else{$ss="$ss ";}
         $j++;
    }
     if($p>$maxp &&  length($longseq)>100){ $good=">$contig\t$loc\t$dir\t$p\t1\n$find\n$sum1\n$ss\n$longseq\n";  $gss=$ss; $maxp=$p; }
###########    gRNA with two insertion
   $sum2="$a1$t2$b1";
   @temp2=split//,$sum2;
   $j=0;  $ss=""; $p=0;
   while($j<($len+1))
    {      $tt="$ar[$j]$temp2[$j]";
         if($ar[$j] eq $temp2[$j]  &&  $ar[$j] ne ""){$ss="$ss*";  $p++;}
           elsif($tt eq "AG" || $tt eq "CT"){$ss="$ss*";  $p++;}
          else{$ss="$ss ";}
         $j++;
    }
     if($p>$maxp &&  length($longseq)>100){ $good=">$contig\t$loc\t$dir\t$p\t2\n$find\n$sum2\n$ss\n$longseq\n";  $gss=$ss; $maxp=$p; }
##########   Contig with one insertion
   $sum11="$a2$t1$b2";
   @temp3=split//,$sum11;
   $j=0;  $ss=""; $p=0;
   while($j<($len+1))
    {  $tt="$temp3[$j]$br[$j]";
        if($br[$j] eq $temp3[$j] &&  $ar[$j] ne ""){$ss="$ss*";   $p++;}
          elsif($tt eq "AG" || $tt eq "CT"){$ss="$ss*";  $p++;}
          else{$ss="$ss ";}
         $j++;
    }
     if($p>$maxp &&  length($longseq)>100){ $good=">$contig\t$loc\t$dir\t$p\t3\n$sum11\n$gRNA\n$ss\n$longseq\n";   $gss=$ss;  $maxp=$p; }
##########   Contig with two insertion
   $sum22="$a2$t2$b2";
   @temp4=split//,$sum22;
   $j=0;  $ss=""; $p=0;
   while($j<($len+1))
    {    $tt="$temp4[$j]$br[$j]";
        if($br[$j] eq $temp4[$j] &&  $ar[$j] ne ""  ){$ss="$ss*";  $p++;}
         elsif($tt eq "AG" || $tt eq "CT"){$ss="$ss*";  $p++;}
          else{$ss="$ss ";}
         $j++;
    }
     if($p>$maxp &&  length($longseq)>100){ $good=">$contig\t$loc\t$dir\t$p\t4\n$sum22\n$gRNA\n$ss\n$longseq\n";   $gss=$ss; $maxp=$p; }
##########
      $i++;
}    

    $gss1=$gss;  $gss2=$gss;  $gss3=$gss; $gss1=~s/\*\*\*\*\*//;  $gss2=~s/\*\*\*\*//;  $gss2=~s/\*\*\*\*//;  $gss3=~s/\*\*\*\*\*\*\*\*//;
     if(length($gss)>(length($gss1)+7) || length($gss)>(length($gss2)+7))
      { if($maxp<($number-3) && ($maxp>($numer/2) || length($gss)>length($gss3) ) ){return $good;} }
}

