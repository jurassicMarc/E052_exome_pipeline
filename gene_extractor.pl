our $gene_list = "/Users/marc/Documents/work/Projects/E052/E052_gene_list_line-endings.csv";

our $vep_dir = "/Users/marc/Documents/work/Projects/E052/variant_annotation";

our (%genes, %variants);

opendir VEP, "$vep_dir" or die $!;

our @vep_files = map {$vep_dir."/".$_} grep {$vep_dir."/".$_ && /.+repeat.+vep/} readdir VEP;

closedir VEP;

open GENES, "$gene_list" or die $!;

while (my $line = <GENES>){
	chomp $line;
	my @genes = split(/\,/, $line);

	my $gene = $line;

	unless (exists $genes{$gene}){
		$genes{$gene} = $gene;
		print "Gene: $genes{$gene}\n";
	}
}

close GENES;

print keys %genes;
print "\n\n";
for my $gene (sort {"\L$a" cmp "\L$b"} keys %genes){
	print "Key: $gene \t Value: $genes{$gene}\n"
}
print "\n";

foreach my $vep (sort {"\L$a" cmp "\L$b"} @vep_files){
	print "VEP_file: $vep\n";
	
	next unless $vep =~ /010/;

	open VEP, "$vep" or die "vep file ($vep) not found, exiting\n";

	my $subject = $vep;
	$subject =~ s/^.+\///g;
	$subject =~ s/\_.+$//;

	my $subject_file = "/Users/marc/Documents/work/Projects/E052/variant_annotation/$subject\_filtered_variants.vep";

	while (my $line = <VEP>){
		chomp $line;
	
		if ($line =~ /^#/){
			next;
			#$line =~ s/^\#//;
			#$variants{$subject}{"header"} .= "\#Gene\t$line";
		}
		
		else{
			if ($line =~ /SYMBOL\=[^\;]+\;/){
				my $gene = $&;
				
				$gene =~ s/SYMBOL\=//;
				$gene =~ s/\;$//;
				
				#if ($gene eq "SMAD5"){
					#print "VEP gene edited: $gene\n";
					#print "gene: $gene\tgenes{gene}: $genes{$gene}\n";
					#print keys %genes; 	
					
					if (exists $genes{$gene}){
						push @{$variants{$subject}{"variants"}{$gene}}, $line;
						#print "Found a variant in GOI ($gene): $line\n";
					}
				#}
				
			}
		}
	};

	open SUBJECT, ">$subject_file" or die "subject file ($subject_file) cannot be opened\n";

	print SUBJECT "$variants{$subject}{header}\n";

	for my $gene (sort {"\L$a" cmp "\L$b"} keys %{$variants{$subject}{"variants"}}){
		foreach my $variant (sort {"\L$a" cmp "\L$b"} @{$variants{$subject}{"variants"}{$gene}}){
			print SUBJECT "$gene\t$variant\n";
		}
	};

	close SUBJECT;

}
