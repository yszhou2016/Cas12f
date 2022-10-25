#Example: perl 3.Indel_Calculate.pl  Sample.R1.fastq.gz   spacer_sequence   reference_sequence

$filein=$ARGV[0];  
$fileout="$filein.output.txt";
$spacer=$ARGV[1];  
$reference=$ARGV[2];

open(FF,">ref.fa"); print FF ">ref\n$reference\n"; close FF;
#`bowtie2-build  ref.fa  ref.fa`;
#`bowtie2  --very-sensitive  -x   ref.fa   -p  3  -U  $filein   --no-unal   -S  $filein.sam`;
`bwa index ref.fa  ref.fa`;
`bwa mem  -A2 -O3 -E1 -t 3 ref.fa  $filein > $filein.sam  `;

$spacer=~tr/atcg/ATCG/;  $len=length($spacer);  $bad_input=0;
$reference=~tr/atcg/ATCG/;  $r=$reference;  $r=~s/$spacer/\t/;
if($reference ne $r){@b=split/\t/,$r; $start=length($b[0]); $end=$start+$len;}
else{ $spacer=~tr/ATCG/TAGC/;  $spacer=reverse($spacer); 
      $r=$reference;    $r=~s/$spacer/\t/;
	  if($reference ne $r){@b=split/\t/,$r; $start=length($b[0]);  $end=$start+$len;}
	  else{$bad_input=1;}
}


open(BB,"$filein.sam");   
open(C1,">temp1");
open(C2,">temp2");
$Edit=0;  $Match=0;  $Deletion=0;  $Insertion=0;  $to=0;
while($line=<BB>)
{     chomp $line;
      @bb=split/\s+/,$line;
      $ll=$bb[5];   $ll2=$bb[5];    $D=0;   $I=0;  $D2=0;  $I2=0; $overlap=0;  
      $ll=~s/M/\t/g;  $ll=~s/D/\t/g; $ll=~s/I/\t/g;
      $ll2=~s/M/M\t/g;  $ll2=~s/D/D\t/g; $ll2=~s/I/I\t/g;
      @ar=split/\s+/,$ll;    @ar2=split/\s+/,$ll2;
      $lenth=@ar;  
      if($bb[3]>1 && $bb[0]!=~/^@/){$to++;   
      if($lenth>1){
          $seq0="_"x($bb[3]-1);  
          $i=0;  $ss=0;  
             $sum="$seq0";
          while($i<$lenth){   
               if($ar2[$i]=~/M$/){   $seq[$i]=substr($bb[9],$ss,$ar[$i]);    $ss=$ss+$ar[$i];  } 
               elsif($ar2[$i]=~/D$/){   $seq[$i]="-"x$ar[$i];  if(($ss+$bb[3]-1)>($end+10) || ($ss+$ar[$i]+$bb[3]-1)<$start){}else{$overlap=1; $D2=1;   $D=$D+$ar[$i];}  }    
               elsif($ar2[$i]=~/I$/){   $seq[$i]=substr($bb[9],$ss,$ar[$i]);  $seq[$i]=~tr/ATGC/atgc/;  
           if(($ss+$bb[3]-1)>($end+10) || ($ss+$ar[$i]+$bb[3]-1)<$start){}else{$overlap=1; $I2=1; $I=$I+$ar[$i]; }    $ss=$ss+$ar[$i]; }    
             $sum="$sum$seq[$i]";
           $i++;  
           }
             if($overlap>0){  $Edit++;  
                        if($D2>0){$Deletion++;}
                        if($I2>0){$Insertion++;}
                 $string="$sum\tDel_len:$D\tIns_len:$I";    #$string=substr($string,$start-20,300);
                        print  C1  $string,"\n";
                 if(not exists $v{$string}){  $v{$string}=1;  }else{  $v{$string}++;  }
                      }
      }else{     
                                     if($ar2[0]=~/M$/ &&   $ar[0]>100){
                                        $seq0="_"x($bb[3]-1);
                                          $seq1=substr($bb[9],0,$ar[0]);      
                                                    $string="$seq0$seq1\tDel_len:0\tIns_len:0";    
                           if(not exists $v{$string}){ $v{$string}=1; }else{ $v{$string}++; }   
                                        print  C2  $string,"\n";    
                                         $Match++;   
                                     }
               }
           }
}
$total=$Edit+$Match;
$RateE=int(10000*$Edit/($total+1))/100;
$RateW=int(100*(100-$RateE))/100;
close BB; close C1;  close C2;

open(AA,"temp1");   open(BB,">temp11");  while($line=<AA>){chomp $line;  @bb=split/\s+/,$line;     $rate=int(10000*$v{$line}/$total)/100;  if(not exists($vv{$line})){print BB $bb[0],"\t","Rate: $rate (%)\t","Count:$v{$line}","\t$bb[1]\t$bb[2]\n";  $vv{$line}=1;  }  }  close AA; close BB;
open(AA,"temp2");   open(BB,">temp22");  while($line=<AA>){chomp $line;  @bb=split/\s+/,$line;       $rate=int(10000*$v{$line}/$total)/100;   if(not exists($vv{$line})){print BB $bb[0],"\t","Rate: $rate (%)\t","Count:$v{$line}","\t$bb[1]\t$bb[2]\n";  $vv{$line}=1;  }  }  close AA; close BB;
`sort  -dk3,3nr  temp11  >  temp11.sorted`;
`sort  -dk3,3nr  temp22  >  temp22.sorted`;
open(AA,"temp22.sorted");  $line=<AA>;   @bb=split/\s+/,$line;  $WT=$bb[0];  close AA;

open(FF,">$fileout");
     print FF "Editing_rate: $RateE %\t","Editing_Reads: $Edit\t","Total_Reads: $total","\t","Deletion_Reads: $Deletion","\t","Insertion_Reads: $Insertion","\n\n>reference\n$reference\n";
    if($bad_input==0){print FF "\ "x($start-8),"Target: ","*"x$len,"\n";}
    else{print FF "Wrong input target sequences!";}
    print FF "$WT\tRate: $RateW (%)\t","Count:$Match      Del_len:0       Ins_len:0\n";  
close FF;
`cat temp11.sorted >> $fileout`;
`rm temp*`;
`rm $filein.sam`;
