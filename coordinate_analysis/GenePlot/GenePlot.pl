use strict;
use warnings;
use Getopt::Std;

my %options;
getopts('g:i:lo:p', \%options);
chomp($options{g},$options{i});
my $GFF = $options{g};
my $queryfile = $options{i};
my $outputpath = $options{o};

my $extension = 'png';
if(defined($options{p})) { $extension = 'pdf'; }

my %start;
my %stop;

open my $in, '<', $GFF;
my $header = <$in>;

while (<$in>) {

	my @columns = split("\t", $_);
	if ($columns[2] eq 'gene') {

		my($gene) = $columns[8] =~ m/ID=(.*?);/;
		$start{$gene} = $columns[3];
		$stop{$gene} = $columns[4];

	}

}

close $in;

open my $query, '<', $queryfile;

my $output_name;
my $gene1;
my $gene2;
my %unique;

while (<$query>) {

chomp $_;
$output_name = $_;
$gene1 = (split('-', $_))[0];
$gene2 = (split('-', $_))[1];

my @coordset = ($start{$gene1},$stop{$gene1},$start{$gene2},$stop{$gene2});
my @sorted_coords = sort { $a <=> $b } @coordset;

my $Start = shift(@sorted_coords);
my $Stop = pop(@sorted_coords);
my $Length = $Stop - $Start;

open my $py, '>', "$outputpath/$output_name.py";
print $py "import matplotlib.pyplot as plt\nfrom dna_features_viewer import GraphicFeature, GraphicRecord\n\n";
#print $py "plt.rcParams[\"axes.linewidth\"] = 1.5\nplt.rcParams[\"axes.linewidth\"] = 1.5\nplt.rcParams[\"font.family\"] = 'Arial'\nplt.rcParams[\"font.weight\"] = 'bold'\nplt.rcParams[\"font.size\"] = 12\n\n";
print $py "plt.rcParams[\"axes.linewidth\"] = 1.5\nfont12 = {'family':'Arial','color':'black','weight':'bold'}\n\n";
#print $py "plt.rcParams[\"axes.linewidth\"] = 1.5\nfont12 = {'family':'Arial','color':'black','weight':'bold','size':12}\n\n";
print $py "features = [\n";

my @annotations;
%unique = ();

open $in, '<', $GFF;
my $header = <$in>;

my $flag = 0;
while (<$in>) {

	chomp $_;
	my @columns = split("\t", $_);
	if ($columns[2] eq 'gene') {
	my($gene) = $columns[8] =~ m/ID=(.*?);/;

		if ($gene eq $gene2) {

			my ($array, $hash) = annotate($_);
			push @annotations, @$array;
			%unique = (%unique, %$hash);
			last;

		}

		if ($flag == 1) {

			my ($array, $hash) = annotate($_);
			push @annotations, @$array;
			%unique = (%unique, %$hash);

		}

		if ($gene eq $gene1) {

			$flag = 1;
			my ($array, $hash) = annotate($_);
			push @annotations, @$array;
			%unique = (%unique, %$hash);

		}

	}

}

close $in;

for my $entry (@annotations) {

	print $py $entry;

}

my $buffer_amt = $Length/100;
my $buffernear = $Start-$buffer_amt;
my $bufferfar = $Length+(2*$buffer_amt);

print $py "]\n\nrecord = GraphicRecord(sequence_length=$bufferfar,first_index=$buffernear, features=features)\nax, _ = record.plot(with_ruler=False)\n";
print $py "ax.figure.savefig('$output_name.$extension', dpi=400, bbox_inches='tight')";
close $py;

chmod 0770, "$outputpath/$output_name.py";
system("python3 $outputpath/$output_name.py");

}

close $query;

sub annotate {

	my $sandwichtop = "\tGraphicFeature(";
	my $sandwichbottom = "),\n";

	my %localunique = %unique;
	my @tannotations;
	my $entry = $_[0];
	my @lines = split("\n", $entry);

	for my $line (@lines) {

		my @columns = split("\t", $line);
		my $type = $columns[2];
		my $beginning = $columns[3];
		my $end = $columns[4];
		my $strand = $columns[6];
		my $name;

		$name = $1;

		if ($type =~ /gene/) {

			my $ID = "$beginning"."$end"."$type";

			unless((exists($unique{$ID}) || (exists($localunique{$ID})))) {

				$localunique{$ID} = 1;

				if ($type =~ /gene/) {

					chomp $name;
					my $string;

					if (defined($options{l})) { $string = "$sandwichtop"."start=$beginning, end = $end, strand=$strand"."1, color=\"#000000\", label=\"$name\", linewidth = 1.5, fontdict = font12"."$sandwichbottom"; }
					else { $string = "$sandwichtop"."start=$beginning, end = $end, strand=$strand"."1, color=\"#000000\", linewidth = 1.5, fontdict = font12"."$sandwichbottom"; }
					push @tannotations, $string;

				}

			}

		}

	}

	return \@tannotations, \%localunique;

}
