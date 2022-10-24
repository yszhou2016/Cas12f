#Example: perl 3.Indel_Caculate.pl  Sample.R1.fastq.gz   spacer_sequence   reference_sequence

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
open(C1,">temp");

print  C1 $reference,"\n";
if($bad_input==0){print  C1 "\ "x($start-8),"Target: ","*"x$len,"\n";}
 else{print C1 "Wrong input target sequences!";}

$Edit=0;  $total=0;  $Deletion=0;  $Insertion=0;
while($line=<BB>)
{     chomp $line;
      @bb=split/\s+/,$line;
      $ll=$bb[5];   $ll2=$bb[5];    $D=0;   $I=0;  $D2=0;  $I2=0; $overlap=0;
      $ll=~s/M/\t/g;  $ll=~s/D/\t/g; $ll=~s/I/\t/g;
      $ll2=~s/M/M\t/g;  $ll2=~s/D/D\t/g; $ll2=~s/I/I\t/g;
      @ar=split/\s+/,$ll;    @ar2=split/\s+/,$ll2;
      $lenth=@ar;  
      if($bb[3]>0 && $bb[0]!=~/^@/){$total++;}
      if($lenth>1){
          $seq0="_"x($bb[3]-1);  
          $i=0;  $ss="";  
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
          print  C1  "$sum\tDel_len:$D\tIns_len:$I","\n";
         }
      }
}
$rate=int(10000*$Edit/($total+1))/100;
close BB; close C1;
open(FF,">$fileout");
print FF "Editing_rate: $rate %\t","Editing_Reads: $Edit\t","Total_Reads: $total","\t","Deletion_Reads: $Deletion","\t","Insertion_Reads: $Insertion","\n\n>reference\n";
close FF;
open(FF,">$fileout.result");
print FF "$filein\tEditing_rate: $rate %\t","Editing_Reads: $Edit\t","Total_Reads: $total","\t","Deletion_Reads: $Deletion","\t","Insertion_Reads: $Insertion","\n\n>reference\n";
close FF;
open(FF,">$fileout.result");
print FF "$filein\tEditing_rate: $rate %\t","Editing_Reads: $Edit\t","Total_Reads: $total","\t","Deletion_Reads: $Deletion","\t","Insertion_Reads: $Insertion","\n";
close FF;
`cat temp >> $fileout`;
`rm temp`;
`rm $file.sam`;
