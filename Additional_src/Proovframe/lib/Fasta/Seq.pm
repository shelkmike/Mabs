package Fasta::Seq;
use warnings;
use strict;

# preference libs in same folder over @INC
use lib '../';

use overload
    'bool' => sub{1},
    '.' => \&cat,
    '""' => \&string;



our $VERSION = '1.0.0';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;

##------------------------------------------------------------------------##

=head1 NAME

Fasta::Seq.pm

=head1 DESCRIPTION

Class for handling FASTA sequences.

=head1 SYNOPSIS

  my $fa = Fasta::Seq->new(">id\natgc);

=back

=cut



##------------------------------------------------------------------------##

=head1 Class Attributes

=cut

our $Base_content_scans = {
	'N' => sub{	return $_[0] =~ tr/'N'// },
	'A' => sub{	return $_[0] =~ tr/'A'// },
	'T' => sub{	return $_[0] =~ tr/'T'// },
	'G' => sub{	return $_[0] =~ tr/'G'// },
	'C' => sub{	return $_[0] =~ tr/'C'// },
};

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

=head2 Add_base_content_scan

=cut

sub Add_base_content_scan{
	my ($class, $patt) = @_;
	return unless defined $patt;
	$Base_content_scans->{$patt} = eval 'sub{ return $_[0] =~ tr/'.$patt.'//; }';
}

=head2 Complement

=cut

sub Complement{
	my ($class, $seq) = @_;
	$seq =~ tr/ATGCatgc/TACGtacg/;
	return $seq;
}

=head2 Reverse_complement

=cut

sub Reverse_complement{
	my ($class, $seq) = @_;
	$seq =~ tr/ATGCatgc/TACGtacg/;
	return scalar reverse $seq;
}

##------------------------------------------------------------------------##


=head1 Constructor METHOD

=head2 new

Create a new FASTA seq object. Either provide a single STRING
 containing one FASTA record or a key => value construct, for which
 either C<seq_head> or C<id> is required.

  my $seq = ">seq1 blub\nATGC\n";
  Fasta::Seq->new($seq);
  # or
  Fasta::Seq->new(
  	id => 'seq1',
  	desc => 'blub',
  	seq => "ATGC",
  );
  # or
  Fasta::Seq->new(
  	seq_head => '>seq1 blub'
  	seq => "ATGC"
  );


=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	# object method -> clone + overwrite
	if($class = ref $proto){
		$self = bless ({%$proto, @_}, $class)
	}else{ # init emtpy obj
		$class = $proto;
		$self = {
				seq_head => '',
				id => '',
				desc => '',
				seq => '',
				byte_offset => undef,
                                window_offset => 0,
			};
	}

	if(@_){
		if(@_%2){ # input is string to split
			my %self;
			my $head;
			($head, $self->{seq}) = split("\n", shift, 2);
			($self->{id}) = $head =~ />?(\S*)/;
			($self{desc}) = $head =~ /[^\S\n](.+)$/;
			$self = {
				%$self,
				%self,
				@_	# overwrite defaults
			};
			$self->{seq_head} = '>'.$self->{id};
			if($self->{desc}){
				$self->{seq_head}.= ' '.$self->{desc}
			}else{
				$self->{desc} = '';
			}

		}else{
			$self = {
				%$self,
				@_	# overwrite defaults
			};

			# make sure, id has not leading '>'
			$self->{id} =~ s/^>// if $self->{id};
			# make sure there is a seq_head entry
			if(!$self->{seq_head} && $self->{id}){
				$self->{seq_head} = '>'.$self->{id};
				$self->{seq_head}.= ' '.$self->{desc} if $self->{desc};
			}elsif($self->{seq_head}){
				my($id,$desc) = $self->{seq_head} =~ m/
					(?:>?(\S*))			# id, >? for records
					(?:[^\S\n]([^\n]+))?		# desc, optional
				/xs;
				$self->{id} = $id;
				$self->{desc} = $desc || '';

			}
		}
		# make sure, head has leading '>'
		substr($self->{seq_head},0,0,'>') unless $self->{seq_head} =~ /^>/;
		chomp($self->{seq});
		chomp($self->{seq_head});
		$self->{seq} =~ tr/\n//d; 	# remove all newlines from seq
	}

#	unless($self->{id} || $self->{seq}){
#		warn "Creation of incomplete FASTA entry: ".
#			($self->{id} ? "sequence missing" : "id missing");
#	}

	return bless $self, $class;
}







##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2 reverse_complement

Reverse complement the sequence.

=cut

sub reverse_complement{
	my $self = shift;
	$self->{seq} = ref($self)->Reverse_complement($self->{seq});
	return $self;
}

=head2 reverse_complement

Complement the sequence.

=cut

sub complement{
	my $self = shift;
	$self->{seq} = ref($self)->Complement($self->{seq});
	return $self;
}


=head2 cat

Concatenate Fasta::Seq object with either another Fasta::Seq object or a
 plain STRING (sequence only). Can be used as Class method as well as Object
methods.
Returns a new object. Keeps the id and other attributes from the first
 provided object.
If the third parameter in Class syntax, the second one in object syntax is
 set to TRUE, the operands are swapped.

Fasta::Seq overloads "." with this method.

  # Class method
  $fab = Fasta::Seq->cat($fa, $fb);
  $fba = Fasta::Seq->cat($fa, $fb, 1); # swap order
  $faATGC = Fasta::Seq->cat($fa, 'ATGC'); # append plain sequence
  $fATGCa = Fasta::Seq->cat('ATGC', $fa); # prepend plain sequence

  # Object method
  $fab = $fa->cat($fb);
  $fba = $fa->cat($fb, 1); # swap order
  $fba = $fa->cat('ATGC', 1); # append + swap -> prepend plain sequence

  # Overload
  $fATGCa = 'ATGC'.$fb; # prepend
  $fa.= $fb; # append by $fb and overwrite $fa

=cut

sub cat{
	my $class = shift unless ref $_[0]; # class usage
	my ($s1, $s2, $swap) = @_;
	unless(ref $s1){
		die 'At least one operand has to be as Fasta::Seq object' unless ref $s2;
		($s1,$s2) = ($s2,$s1);
		$swap = $swap ? 0 : 1; # toggle swap
	}

	my $re = $s1->new; # clone
	if($swap){
		$re->seq(ref $s2
			? $s2->seq.$re->seq
			: $s2.$re->seq
		);
	}else{
		$re->seq(ref $s2
			? $re->seq.$s2->seq
			: $re->seq.$s2
		);
	}

	return $re;
}


=head2 base_content

=cut

sub base_content{
	my $self = shift;
	my $patt = shift;
	__PACKAGE__->Add_base_content_scan($patt) unless exists $Base_content_scans->{$patt};
	return &{$Base_content_scans->{$patt}}($self->seq)
}


=head2 substr_seq

Substr sequence. Takes the same parameter as perls C<substr>, either plain
 or as a LIST of ARRAYREFS, each containing a set of parameter to allow for
 multiple operations at once.

In the first case the objects sequence is modified, the description appended
 by SUBSTR:<OFFSET>,<LENGTH> and the modified object is returned.

In the second case, the object itself is not modified, but a LIST of modified
 clones will be returned. The ids are appended by <.CLONECOUNTER>.
Methods like qual_low or qual_lcs provide these kind of LIST of ARRAYREFS.

  $fq->substr_seq(5);   	# nt 5 to end
  $fq->substr_seq(5,-5);    # nt 5 to fifth nt from the end
  $fq->substr_seq(-5)		# last 5 nts
  $fq->substr_seq(-100,50)  # 50 nts, starting 100 from end
  $fq->substr_seq(10,5,"AAAA") # replace 5 nts starting at pos 10
    with 4 "A"s and replace corresponding qual values.
  ($fq1, $fq2) = $fq->substr_seq([5,100], [200,500]);
    # nts 5 to 100 and 200 to 500


=cut

sub substr_seq{
	my $self = shift;

	if(! ref $_[0] || @_ == 1){

		my $fa = $self->new; # clone

		my ($o, $l, $r) = ref $_[0] ? @{$_[0]} : @_;

		die __PACKAGE__."::substr_seq: Not enougth arguments\n"
			unless defined ($o);

		# replace
		if(defined $r){
			$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
			substr($fa->{seq}, $o, $l, $r);
                    }elsif(defined $l){
			$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
			$fa->seq( substr($fa->{seq}, $o, $l) );
		}else{
			$fa->desc_append(sprintf("SUBSTR:%d", $o));
			$fa->seq( substr($fa->{seq}, $o) );
		}
		return $fa;
	}else{
		my @new_fas;
		my $clone_c = 0;
		foreach (@_){
			$clone_c++;
			my $fa = $self->new; # clone
			$fa->id($fa->id.".$clone_c");

			my ($o, $l, $r) = @$_;

			die __PACKAGE__."::substr_seq: Not enougth arguments\n"
				unless defined ($o);

			# replace
			if(defined $r){
				$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
				substr($fa->{seq}, $o, $l, $r);
			}elsif(defined $l){
				$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
				$fa->seq( substr($fa->{seq}, $o, $l) );
			}else{
				$fa->desc_append(sprintf("SUBSTR:%d", $o));
				$fa->seq( substr($fa->{seq}, $o) );
			}
			push @new_fas, $fa;
		}
		return @new_fas;
	}

}

=head2 desc_append

Append description by given attributes. Attributed are appended with
 whitespace delimiter.

=cut

sub desc_append{
	my ($self, @attr) = @_;
	my $cur_desc = $self->desc();
	$self->desc(
		$cur_desc
			?	join(" ", $cur_desc, @attr)
			:	join(" ", @attr)
	);
	return $self->{desc};
}

=head2 string([LINEWIDTH])

Get entire sequence as FASTA string. Takes line_width => INT, default 80, 0 for single line, and fastq => CHAR to return FASTQ string.

=cut

sub string{
    my ($self, $lw) = (@_);
    if (! defined($self->{seq}) || $self->{seq} eq ""){
        warn __PACKAGE__."::string: >", $self->id, " has empty sequence. Returning empty sting";
        return "";
    }
    if ( @_ != 2 && defined $lw) { # not single/no option: overloaded (undef) or single opt line width
        my %p = (line_width => 80, @_[1..$#_]);
        if ($p{fastq}) { # return FASTQ string
            return '@'.(substr ($self->{seq_head}, 1))."\n".
                $self->{seq}."\n+\n".
                ($p{fastq} x length($self->seq))."\n";
        }else {
            $lw = $p{line_width};
        }
    }

    $lw //= 80; # default line width is 80
    if($lw){
        my $s = "";
        $s.= $_."\n" for unpack "(A$lw)*", $self->{seq};
        return $self->{seq_head}."\n".$s;
    }else{
        return $self->{seq_head}."\n".$self->{seq}."\n";
    }
}

=head2 next_window(SIZE, [SHIFT=SIZE+1])

Move along sequence, returning subsequences of given SIZE. Last window will be >= SIZE, unless the entire sequence is < SIZE;

=cut

sub next_window{
    my ($self, $size, $shift, $overhang) = @_;
    $shift//= $size;
    $overhang//= "keep";
    die "size and shift need to be > 0" unless ($size > 0 and $shift > 0);

    if(!defined($self->{window_offset})){
        $self->{window_offset} = 0;
        return undef
    };

    my $offset = $self->{window_offset};
    if($offset+$shift+$size > $self->length){ # next win is overhang
        if ($overhang eq "keep") { # keep slicing overhang
            $self->{window_offset} = $offset+$shift;
            $self->{window_offset} = undef if $self->{window_offset} >= $self->length;
            $self->substr_seq($offset, $size);
        }elsif ($overhang eq "drop") { # return window, ignore overhang
            $self->{window_offset} = undef;
            $self->substr_seq($offset, $size);
        }elsif ($overhang eq "merge") { # return window + overhang
            $self->{window_offset} = undef;
            $self->substr_seq($offset);
        }else {
            die "unknown mode $overhang for handling overhang window";
        }
    }else {
        $self->{window_offset} =  $offset+$shift;
        $self->substr_seq($offset, $size);
    }
}

##------------------------------------------------------------------------##

=head1 Accessor METHODS

=head2 seq_head

Get/Set the seq_head. Also updates id and desc.

=cut


sub seq_head{
	my ($self, $seq_head) = @_;
	if ($seq_head){
		$self->{seq_head} = $seq_head;
		# make sure, head has leading '>'
		substr($self->{seq_head},0,0,'>') unless $self->{seq_head} =~ /^>/;
		# reset id cache if head is changed
		my($id,$desc) = shift =~ m/
			(?:>?(\S+))			# id, >? for records
			(?:\s([^\n]+))?		# desc, optional
		/xs;
		$self->{id} = $id;
		$self->{desc} = $desc;
	};
	return $self->{seq_head};
}

=head2 seq

Get/Set the seq.

=cut

sub seq{
	my ($self, $seq) = @_;
	if(defined $seq){
                $seq =~tr/\n//d;
                warn __PACKAGE__."::seq: >", $self->id, " set to zero-length sequence. Not a valid Fasta record\n" if $seq eq "";
		$self->{seq} = $seq
	};
	return $self->{seq};
}

=head2 byte_offset

Get/Set the byte_offset.

=cut

sub byte_offset{
	my ($self, $byte_offset) = @_;
	$self->{byte_offset} = $byte_offset if defined($byte_offset);
	return $self->{byte_offset};
}


=head2 id

Get/set the seqs id. Also updates seq_head.

=cut

sub id{
	my ($self, $id) = @_;
	if(defined $id){
		$self->{id} = $id;
		# reset seq_head
		# make sure, desc is cached
		$self->{seq_head} = $self->{desc}
			? '>'.$self->{id}.' '.$self->{desc}
			: '>'.$self->{id};
	}
	return $self->{id};
}

=head2 desc

Get/set the seqs description, also updates seq_head.

=cut

sub desc{
	my ($self, $desc) = @_;
	if(defined $desc){
		$self->{desc} = $desc;
		# reset seq_head
		# make sure, desc is cached
		$self->{seq_head} = $self->{desc}
			? '>'.$self->{id}.' '.$self->{desc}
			: '>'.$self->{id}; # desc is '' ->remove desc
	}
	return $self->{desc};
}

=head2 length

Get seq length.

=cut

sub length{
    length($_[0]->{seq});
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
