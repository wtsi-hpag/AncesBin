#!/usr/bin/perl -w

$rdname1 = "";
$rdname2 = "";
$contig1 = "";
$contig2 = "";
$dirtag1 = 1;
$dirtag2 = 1;
$mscore1 = 1;
$mscore2 = 1;
$n_Hap1 = 0;
$n_Hap2 = 0;
$m_Hap1 = 0;
$m_Hap2 = 0;

while ($line = <>) 
{
     chomp ($line); 
     @fields = split(/\s/,$line);
     if($rdname1 eq $fields[0])
     {
       $rdname2 = $fields[0];
       $dirtag2 = $fields[1];
       $contig2 = $fields[2];
       $mscore2 = $fields[3];
       $n_Hap1 = 0;
       $n_Hap2 = 0;
       $m_Hap1 = 0;
       $m_Hap2 = 0;
       if(substr($contig1,0,2) =~ /MA/)
       {
         $n_Hap1 = $n_Hap1+100;
         $m_Hap1 = $m_Hap1+$mscore1;
        
#           printf "%s %d %s %s || %s %d %s %s\n",$rdname1,$dirtag1,$contig1,$mscore1,$rdname2,$dirtag2,$contig2,$mscore2;
       }
       else
       {
         $n_Hap2 = $n_Hap2+100;
         $m_Hap2 = $m_Hap2+$mscore1;
       }

       if(substr($contig2,0,2) =~ /MA/)
       {
         $n_Hap1 = $n_Hap1+100;
         $m_Hap1 = $m_Hap1+$mscore2;
        
#           printf "%s %d %s %s || %s %d %s %s\n",$rdname1,$dirtag1,$contig1,$mscore1,$rdname2,$dirtag2,$contig2,$mscore2;
       }
       else
       {
         $n_Hap2 = $n_Hap2+100;
         $m_Hap2 = $m_Hap2+$mscore2;
       }

       if(($m_Hap2 == 0)&&($m_Hap1 == 0))
       {
          printf("%s MAT\n",$rdname1);
          printf("%s PAT\n",$rdname1);
       }
       elsif(($n_Hap2 == 100)&&($n_Hap1 == 100))
       {
         if($m_Hap1>$m_Hap2)
         {
           printf("%s MAT\n",$rdname1);
         }
         elsif($m_Hap1<$m_Hap2)
         {
           printf("%s PAT\n",$rdname1);
         }
         else
         {
           printf("%s MAT\n",$rdname1);
           printf("%s PAT\n",$rdname1);
         }
       }
       else
       {
         if($n_Hap1>$n_Hap2)
         {
           printf("%s MAT\n",$rdname1);
         }
         elsif($n_Hap2 > 0)
         {
           printf("%s PAT\n",$rdname1);
         }
       }
       $n_contig++;
     }
     else
     {
       $n_contig = 1;
     }
     $rdname1 = $fields[0];
     $dirtag1 = $fields[1];
     $contig1 = $fields[2];
     $mscore1 = $fields[3];
} 


__END__
