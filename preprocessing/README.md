# Preprocessing Scripts for Cross-Species Benchmarks

This directory contains species-specific preprocessing scripts to generate benchmark datasets.

## Species-Specific Preprocessing

### Human (SpliceVarDB)
**Location**: `human/`

**Required input**:
- Reference FASTA: Ensembl r112 (GRCh38)
  - Source: https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/
- GTF Annotation: GENCODE Version 38
  - Source: GENCODE
- sQTL/Variant Data: SpliceVarDB
  - Source: https://compbio.ccia.org.au/splicevardb/

**Output**: 23,877 SNVs (12,602 splice-altering, 11,275 normal)

---

### Rat (RatGTEx)
**Location**: `rat/`

**Required input**:
- Reference FASTA: Ensembl r84 (Rnor_6.0)
  - Source: https://ftp.ensembl.org/pub/release-84/fasta/rattus_norvegicus/dna/
- GTF Annotation: Ensembl r84
  - Source: https://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/
- sQTL Data: RatGTEx v1
  - Source: https://ratgtex.org/download/v1/
- Ensembl VCF: Ensembl r84 (4,709,172 variants)
  - Source: https://ftp.ensembl.org/pub/release-84/variation/vcf/rattus_norvegicus/

**Output**: 28,120 SNVs (14,060 splice-altering, 14,060 normal)

---

### Pig (PigGTEx)
**Location**: `pig/`

**Required input**:
- Reference FASTA: Ensembl r102 (Sscrofa11.1)
  - Source: https://ftp.ensembl.org/pub/release-102/fasta/sus_scrofa/dna/
- GTF Annotation: Ensembl r102
  - Source: https://ftp.ensembl.org/pub/release-102/gtf/sus_scrofa/
- sQTL Data: PigGTEx (FarmGTEx)
  - Source: https://www.farmgtex.org/
- Ensembl VCF: Ensembl r102 (3,847,291 variants)
  - Source: https://ftp.ensembl.org/pub/release-102/variation/vcf/sus_scrofa/

**Output**: 26,358 SNVs (13,179 splice-altering, 13,179 normal)

---

### Chicken (ChickenGTEx)
**Location**: `chicken/`

**Required input**:
- Reference FASTA: Ensembl r102 (GRCg6a)
  - Source: https://ftp.ensembl.org/pub/release-102/fasta/gallus_gallus/dna/
- GTF Annotation: Ensembl r102
  - Source: https://ftp.ensembl.org/pub/release-102/gtf/gallus_gallus/
- sQTL Data: ChickenGTEx portal
  - Source: https://ngdc.cncb.ac.cn/chickengtex/
- Ensembl VCF: Ensembl r102 (~2.5M variants)
  - Source: https://ftp.ensembl.org/pub/release-102/variation/vcf/gallus_gallus/

**Output**: 25,000 SNVs (12,500 splice-altering, 12,500 normal)

