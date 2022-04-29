$file=$ARGV[0];  # genome fasta
$contig=$ARGV[1];  # Contig name of Cas protein
$start=$ARGV[2];   # Start codon position of Cas protein in Contig
$end=$ARGV[3];     # end codon position of Cas protein in Contig
$DR=$ARGV[4];      # direct repeat of Cas protein


######################
open(AA,"$file");
while($line=<AA>)
{  $line2=<AA>;
   chomp $line;  chomp $line2;  $ll=$line;  $ll=~s/>//;
   @bb=split/\s+/,$ll;
   $h{$bb[0]}=$line2;
}
close AA;
######################
 
if($end>$start){ $line=substr($h{$contig},$end,300); }
else{       $line=substr($h{$contig},$end-300,300);
    $line=~tr/ATCGatcg-/TAGCTAGC-/;   $line=reverse($line);
}

##################################################
$number=22;  
$ii=length($DR)-$number;   
$grna=substr($DR,$ii,$number);

$sca=200; #### extend length

$alen=length($line);
$ii=0;
open(BB,">$file.tracrRNA");
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
      print BB $ff1,$ff2; 
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
     if($p>$maxp &&  length($longseq)>70){ $good=">$contig\t$loc\t$dir\t$p\t1\n$find\n$sum1\n$ss\n$longseq\n";  $gss=$ss; $maxp=$p; }
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
     if($p>$maxp &&  length($longseq)>70){ $good=">$contig\t$loc\t$dir\t$p\t2\n$find\n$sum2\n$ss\n$longseq\n";  $gss=$ss; $maxp=$p; }
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
     if($p>$maxp &&  length($longseq)>70){ $good=">$contig\t$loc\t$dir\t$p\t3\n$sum11\n$gRNA\n$ss\n$longseq\n";   $gss=$ss;  $maxp=$p; }
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
     if($p>$maxp &&  length($longseq)>70){ $good=">$contig\t$loc\t$dir\t$p\t4\n$sum22\n$gRNA\n$ss\n$longseq\n";   $gss=$ss; $maxp=$p; }
##########
      $i++;
}    #$gss=~s/ \* /   /g;
     #$gss1=$gss;  $gss2=$gss;  $gss3=$gss; 
     #$gss1=~s/\*\*\*  \*\*\*\*\*//;  $gss1=~s/\*\*\*\*  \*\*\*\*//;  $gss1=~s/\*\*\*\*\*  \*\*\*//;
     #$gss1=~s/\*\*\*\*\*//;  $gss1=~s/\*\*\*\*//;
     #$gss3=~s/\*\*\*\*\*\*\*//;     $gss3=~s/\*\*\*\*\*\*\*\*//;   
     #$gss3=~s/\*\*\*\*\*\*\*\*//;    $gss3=~s/\*\* \*\*\*\*\*//;  
     #$gss3=~s/\*\*\* \*\*\*\*//;    $gss3=~s/\*\*\* \*\*\*\*//;   
     #$gss3=~s/\*\*\*\* \*\*\*//;    $gss3=~s/\*\*\*\*\* \*\*//;   
     #$gss3=~s/\*\*\*\*\*\*\*\*//;
      
     #if($maxp<($number-4) && ($maxp> ($numer/2)+2))
      # { if(length($gss)>(length($gss1)+6) || length($gss)>length($gss2) || length($gss)>(length($gss3)+6))
       #  {return $good;} 
       #}
    $gss1=$gss;  $gss2=$gss;  $gss3=$gss; $gss1=~s/\*\*\*\*\*//;  $gss2=~s/\*\*\*\*//;  $gss2=~s/\*\*\*\*//;  $gss3=~s/\*\*\*\*\*\*\*\*//;
     if(length($gss)>(length($gss1)+7) || length($gss)>(length($gss2)+7))
      { if($maxp<($number-4) && ($maxp>($numer/2) || length($gss)>length($gss3) ) ){return $good;} }

}

