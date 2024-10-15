# Building a Custom SnpEff Database

## Table of Contents
1. Introduction
2. Important Notes: Reading the SnpEff Documentation
3. Prerequisites: Conda Installation
4. Using Pre-built SnpEff Databases
5. When to Build a Custom SnpEff Database
6. Formatting and Errors: Handling GFF Files
7. Verifying Mutations with IGV

---

## 1. Introduction

In this tutorial, we will guide you through the process of building a custom SnpEff database. The steps outlined here are meant for users who need to annotate genomic variants using their own reference genome or custom annotations that may not be available in the pre-built SnpEff databases.

## 2. Important Notes: Reading the SnpEff Documentation

Before proceeding with this tutorial, it is essential that you carefully read the official SnpEff documentation. The official guide provides valuable insights and details on how SnpEff operates, which can help avoid common pitfalls during the database-building process. You can find the official documentation at the following URL:

[SnpEff Documentation](https://pcingola.github.io/SnpEff/snpeff/introduction/)

Taking the time to understand this material will improve the accuracy and reliability of your custom annotations.

## 3. Prerequisites: Conda Installation

This tutorial assumes that you have already installed SnpEff using Conda. If you have not installed it yet, please follow the instructions provided in the official SnpEff documentation or execute the following command to install it:

```bash
conda install -c bioconda snpeff
```

## 4. Using Pre-built SnpEff Databases

If your reference genome is already available in a pre-built SnpEff database, you can simply use the `-e` option followed by the name of the database. This allows you to avoid building a custom database. For example:

```bash
snpeff -e GRCh38.99 input.vcf > annotated.vcf
```

This command will use the pre-built GRCh38.99 database for variant annotation.

## 5. When to Build a Custom SnpEff Database

There are situations where you may need to build your own SnpEff database. These include:

- Your reference genome is not available in the existing SnpEff databases.
- You want to use a custom annotation that you have created or downloaded.

In these cases, this tutorial will guide you through the necessary steps to build your own database from scratch.

## 6. Formatting and Errors: Handling GFF Files

When building a custom SnpEff database, one of the most critical components is ensuring that your GFF (General Feature Format) file is correctly formatted. GFF files are used to describe gene structures and other annotations in your reference genome.

However, GFF files can vary greatly in format, and this can lead to errors during the SnpEff database build process. In some cases, errors may not be immediately obvious, but incorrect formatting can lead to unexpected or incorrect annotation results.

### Common Issues with GFF Files

- Missing or malformed attributes.
- Incorrect gene model structures.
- Formatting discrepancies that are tolerated by some tools but not by SnpEff.

To avoid these issues, it is essential to validate your GFF files using dedicated tools, or to carefully inspect the structure before proceeding with the database build.

## 7. Verifying Mutations with IGV

Even if your SnpEff database is built correctly and no errors are reported, it is highly recommended that you verify the results for high-impact mutations manually. The SnpEff annotation process may flag certain mutations as high-impact, but the actual effect of the mutation on gene function should be validated using genome visualization tools like IGV (Integrative Genomics Viewer).

IGV allows you to visualize the reference genome alongside your variant data, helping you to understand how the variants affect gene structures and functions.

You can download IGV and learn more about it from the following link:

[IGV](https://igv.org)
