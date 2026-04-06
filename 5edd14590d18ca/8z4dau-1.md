---
title: "📰 リガンドの準備"
free: false
---

モデルのリガンドはCMNPD28986（TOP-2)です。

まずはPyMOLを用いて、分子ドッキングしたファイルを開いてください。
※ここでは違う化合物を用いていますが、やり方は同じです。

以下のようにAction→hydrogens→→addを押してください。
![](https://storage.googleapis.com/zenn-user-upload/663e1603fa00-20250924.png)

以下のように水素が付加されます。

![](https://storage.googleapis.com/zenn-user-upload/dd10f131ac96-20250924.png)

左上のFile→Export Molecule...を押し、
以下の画面が出てくると思うので、saveを押してください。

![](https://storage.googleapis.com/zenn-user-upload/7a08751bef41-20250924.png)

その後、`ligand.mol2`としてmol2ファイルで保存してください。

GROMACSではAMBER力場を直接使うことはできないため、リガンドに対してGAFF（General Amber Force Field）を適用するためにACPYPEを使用します。

ligand.mol2の中身をテキストファイルで開いてみてみます。
`# created with PyMOL 3.1.0`
の部分があると、のちにエラーになるので、消去してください。

```
# created with PyMOL 3.1.0　←ここは消去してください。
@<TRIPOS>MOLECULE
CMNPD28986
 29 30 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C           6.9820   26.0100  -21.2750 C.ar    1  UNL1       -0.0198
      2 C           7.5840   25.5080  -20.1250 C.ar    1  UNL1       -0.0577
      3 C           6.9120   25.5620  -18.9070 C.ar    1  UNL1       -0.0540
      4 C           5.6260   26.1050  -18.8120 C.ar    1  UNL1       -0.0164
      5 C           5.0030   26.6070  -19.9900 C.ar    1  UNL1        0.0266
      6 C           5.6980   26.5570  -21.2150 C.ar    1  UNL1        0.1228
      7 C           5.0430   26.1520  -17.4450 C.2     1  UNL1       -0.0556
      8 C           3.8820   26.6900  -17.0380 C.2     1  UNL1       -0.0526
      9 C           2.8460   27.2920  -17.9410 C.3     1  UNL1        0.1031
     10 O           3.4560   28.0730  -18.9240 O.3     1  UNL1       -0.3653
     11 C           3.6140   27.1780  -19.9790 C.3     1  UNL1        0.0775
     12 C           1.8320   28.1280  -17.1560 C.3     1  UNL1        0.0821
     13 C           2.4600   29.3820  -16.5690 C.3     1  UNL1       -0.0362
     14 O           1.2960   27.3660  -16.1180 O.3     1  UNL1       -0.3892
     15 H           1.3020   26.4290  -16.3680 H       1  UNL1        0.2099
     16 O           5.1350   27.0310  -22.3410 O.3     1  UNL1       -0.5064
     17 H           4.2130   26.9390  -22.4770 H       1  UNL1        0.2921
     18 H           7.5044   25.9769  -22.2082 H       1  UNL1        0.0654
     19 H           8.5635   25.0805  -20.1770 H       1  UNL1        0.0619
     20 H           7.3877   25.1813  -18.0274 H       1  UNL1        0.0624
     21 H           5.6311   25.6927  -16.6781 H       1  UNL1        0.0621
     22 H           3.6798   26.6893  -15.9873 H       1  UNL1        0.0599
     23 H           2.3158   26.4890  -18.4089 H       1  UNL1        0.0666
     24 H           3.4375   27.6892  -20.9022 H       1  UNL1        0.0611
     25 H           2.9148   26.3773  -19.8573 H       1  UNL1        0.0611
     26 H           1.0654   28.4249  -17.8409 H       1  UNL1        0.0619
     27 H           3.2307   29.1059  -15.8800 H       1  UNL1        0.0255
     28 H           2.8801   29.9726  -17.3562 H       1  UNL1        0.0255
     29 H           1.7112   29.9501  -16.0576 H       1  UNL1        0.0255
@<TRIPOS>BOND
     1    16    17    1
     2     6    16    1
     3     1     6   ar
     4     1     2   ar
     5     5     6   ar
     6     2     3   ar
     7     5    11    1
     8     4     5   ar
     9    10    11    1
    10     9    10    1
    11     3     4   ar
    12     4     7    1
    13     9    12    1
    14     8     9    1
    15     7     8    2
    16    12    13    1
    17    12    14    1
    18    14    15    1
    19     1    18    1
    20     2    19    1
    21     3    20    1
    22     7    21    1
    23     8    22    1
    24     9    23    1
    25    11    24    1
    26    11    25    1
    27    12    26    1
    28    13    27    1
    29    13    28    1
    30    13    29    1

```

[gromacsのチュートリアル](https://mdtutorials.com/gmx/complex/02_topology.html)によると、この@<TRIPOS>BONDが昇順になっていないとエラーが出るらしいので、以下のスクリプトを  `sort_mol2_bonds.pl`　として作って、実行します。

```perl
#!/usr/bin/perl

use strict;

# sort_mol2_bonds.pl - a script to reorder the listing in a .mol2 @<TRIPOS>BOND
# section so that the following conventions are preserved:
#   1. Atoms on each line are in increasing order (e.g. 1 2 not 2 1)
#   2. The bonds appear in order of ascending atom number
#   3. For bonds involving the same atom in the first position, the bonds appear
#       in order of ascending second atom
#
# Written by: Justin Lemkul (jalemkul@vt.edu)
#
# Distributed under the GPL-3.0 license

unless (scalar(@ARGV)==2)
{
    die "Usage: perl sort_mol2_bonds.pl input.mol2 output.mol2\n";
}

my $input = $ARGV[0];
my $output = $ARGV[1];

open(IN, "<$input") || die "Cannot open $input: $!\n";
my @in = <IN>;
close(IN);

# test for header lines that some scripts produce
unless($in[0] =~ /TRIPOS/)
{
    die "Nonstandard header found: $in[0]. Please delete header lines until the TRIPOS molecule definition.\n";
}

open(OUT, ">$output") || die "Cannot open $output: $!\n";

# get number of atoms and number of bonds from mol2 file
my @tmp = split(" ", $in[2]);
my $natom = $tmp[0];
my $nbond = $tmp[1];

# check
print "Found $natom atoms in the molecule, with $nbond bonds.\n";

# print out everything up until the bond section
my $i=0;
while (!($in[$i] =~ /BOND/))
{
    print OUT $in[$i];
    $i++;
}

# print the bond section header line to output
print OUT $in[$i];
$i++;

# read in the bonds and sort them
my $bondfmt = "%6d%6d%6d%5s\n";
my @tmparray;

# sort the bonds - e.g. the one that has the
# lowest atom number in the first position and then the
# lowest atom number in the second position (swap if necessary)
for (my $j=0; $j<$nbond; $j++)
{
    my @tmp = split(" ", $in[$i+$j]);
    # parse atom numbers
    my $ai = $tmp[1];
    my $aj = $tmp[2];
    # reorder if second atom number < first
    if ($aj < $ai)
    {
        $ai = $tmp[2];
        $aj = $tmp[1];
    }
    # store new lines in a temporary array
    $tmparray[$j] = sprintf($bondfmt, $tmp[0], $ai, $aj, $tmp[3]);
}

# loop over tmparray to find each atom number
my $nbond = 0;
for (my $x=1; $x<=$natom; $x++)
{
    my @bondarray;
    my $ntmp = scalar(@tmparray);
    for (my $b=0; $b<$ntmp; $b++)
    {
        my @tmp = split(" ", $tmparray[$b]);
        if ($tmp[1] == $x)
        {
            push(@bondarray, $tmparray[$b]);
            splice(@tmparray, $b, 1);
            $ntmp--;
            $b--;
        }
    }

    if (scalar(@bondarray) > 0) # some atoms will only appear in $aj, not $ai
    {
        my $nbondarray = scalar(@bondarray);
        if ($nbondarray > 1)
        {
            # loop over all bonds, find the one with lowest $aj
            # and then print it
            for (my $y=0; $y<$nbondarray; $y++)
            {
                my @tmp2 = split(" ", $bondarray[$y]);
                my $tmpatom = $tmp[2];
                my $lowindex = 0;
                if ($tmp2[2] < $tmpatom)
                {
                    $lowindex = $y;
                }
                my $keep = splice(@bondarray, $lowindex, 1);
                $y--;
                $nbondarray--;
                my @sorted = split(" ", $keep);
                $nbond++;
                printf OUT $bondfmt, $nbond, $sorted[1], $sorted[2], $sorted[3];
            }
        }
        else
        {
            $nbond++;
            my @tmp2 = split(" ", $bondarray[0]);
            printf OUT $bondfmt, $nbond, $tmp2[1], $tmp2[2], $tmp2[3];
        }
    }
}

close(OUT);

exit;
```

以下を実行してください。

```perl
 perl sort_mol2_bonds.pl ligand.mol2 ligand_fix.mol2
```

続いて、

```bash
 acpype -i ligand_fix.mol2
```
を実行してください。

すると、`ligand_fix.acpype` というフォルダが出てきており、
その中に様々なファイルがあると思います。

使用していくのは赤枠で囲った.groと.itpファイルです。

![](https://storage.googleapis.com/zenn-user-upload/80340d9c9f48-20250924.png)

これでリガンドの準備は終了です！
