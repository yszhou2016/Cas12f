$file=$ARGV[0];      #   The fasta file of Cas proteins
#############################
# perl  1.Cas12f-Finder.pl   Cas.pep.fa
#############################
$file1=$file;  $file1=~s/\.fasta/\./;  $file1=~s/\.fa\./\./;
open(FASTA,"$file");
open(FA,">$file1.tmp");
$s=0;
while($line=<FASTA>)
{  chomp $line;
    if($line=~/^>/)
   { if($s>0){print FA "\n";}
     print FA $line,"\n";
     $s=1;
   }
    else{$line=~s/\r//g; print FA $line;}
 }
 print FA "\n";
close FASTA; close FA;
######################################
open(AA,"$file1.tmp");
open(BB,">$file1.Cas12f.fa");
#open(CC,">$file1.Cas12f.imfor");
while($line1=<AA>)
{  $line2=<AA>;
   chomp $line1;  chomp $line2;
   $head=substr($line2,1,100);
   if($line2=~/^M/ && length($line2)<700 && length($line2)>350 && not exists($h{$head}))
       {
      @bb=split//,$line2;
     $i=0; $len=@bb;
     $s1=0;  $s2=0;  $s3=0;  $s4=0;  $s5=0;
     #print CC $line1,"\t";  #####
     while($i<$len)
     {  
        if($s1==0){
        if($bb[$i] eq "G" || $bb[$i] eq "S" || $bb[$i] eq "A")
           {  if($bb[$i+2] eq "D")
               {  if($bb[$i+4] eq "G" || $bb[$i+4] eq "N")
                   { $s1=1;   #print CC "RuvC-I","\t";  #####    
                   }
               }
           }
        }
        #################
        if($s1>0 && $s2==0)
         {
           if($bb[$i] eq "V" || $bb[$i] eq "I" || $bb[$i] eq "P" || $bb[$i] eq "L")
           {  if($bb[$i+3] eq "E" )
               {  $s2=1;
                   #print CC "RuvC-II-1","\t";  ##### 
               }
           }
        }
        #################
        if($s2>0 ){
           {  if($bb[$i] eq "T" || $bb[$i] eq "S")
               {  if($bb[$i+1] eq "S"   &&   $bb[$i+4] eq "C" &&   $bb[$i+7] eq "C"  )
                   {  $s3=1;  $i=$i+5;
                     #print CC "RuvC-II-2_Zn_Finger","\t";  ##### 
                   }
               }
           }
        }
        if($s3>0 &&  $s4==0)
             {    if($bb[$i] eq "C"   &&     $bb[$i+3] eq "C"  )
                  {$s4=1;   # print CC "Zn_Finger2","\t"; 
                                   }
             }
        #################
        if($s4>0 && $s5==0){
        if($bb[$i] eq "D" )
           {     if($bb[$i+3] eq "A" || $bb[$i+3] eq "G"  ||   $bb[$i+3] eq "V")
                   {    if($bb[$i+4] eq "A" || $bb[$i+4] eq "S")
                        { $s5=1;
                          #print CC "RuvC-III","\t";  ##### 
                        }
                   } 
           }
         }
        #################   
        $i++;
     }
     print CC "\n";
     if($s4>0  ){print BB $line1,"\n",$line2,"\n";}
    $h{$head}=1;
  }
}
close AA; close BB;  #close CC;

`rm  $file1.tmp`;

