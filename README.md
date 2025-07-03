
File has lines that look like this

```
+ GCGGTTTTTGGTGAGTTTAT..CAACTTTCAGCGGCTCCATC 1130400.0 AGAGCTCCATGTAAGTTAGT..CGAATTTCAGTGCTCCACGC 1778380.0
```

The `+` indicates this is a plus strand gene. Therefore, the splice sites will
generally follow the GT-AG rule. On the negative strand, the sequence is CT-AC.
The splice dontor concensus is `GTAAG`. The splice acceptor concensus is
`TTTCAG`.


The first part before the .. is 10 nt exon followed by 10 nt intron. You can see the start of the intron with the canonical GT.

```
GCGGTTTTTGGTGAGTTTAT
==========----------
  exon       intron
```

The second part after the .. is 10+10 again. The .. represents the middle of
the intron whose length is quite variable.

```
CAACTTTCAGCGGCTCCATC
----------==========
  intron    exon
```

Each line represents adjacent introns. The number after the exon-intron
sequence is the expression level. Here, the intron to the left has a much lower
expression level than the one on the right.

```
+ exon-intron 1130400.0 exon-intron 1778380.0
```
