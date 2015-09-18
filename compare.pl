#!/usr/bin/perl

my $inAlign = 0;

while(<>)
{
    if(/^Alignment\.\.\./)
    {
        $_=<>;
        last;
    }
}

while(<>)
{
    if(/^Score/)
    {
        last;
    }

    my $line1 = $_;
    my $line2 = <>;
    $_ = <>;
    
    Compare($line1, $line2);
}

sub Compare
{
    my($line1, $line2) = @_;

    chomp $line1;
    chomp $line2;

    my @chars1 = split(//, $line1);
    my @chars2 = split(//, $line2);

    my $compare = "";

    for(my $i=0; $i<scalar(@chars1); $i++)
    {
        if($chars1[$i] eq $chars2[$i])
        {
            $compare .= ":";
        }
        else
        {
            $compare .= "X";
        }
    }
    print "$line1\n";
    print "$compare\n";
    print "$line2\n\n";

}
