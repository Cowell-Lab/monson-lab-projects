# MonsonLab MAV Project

Project ID: 16526c92-9ffb-4e1a-96ce-f9ed77baa3a6

This project consists of files from multiple sequencing libraries.

+ Library 4, 6, 8, 9, 10, 23 (43 repertoires)
+ Library 13, 15, 16, 20, (13 repertoires)
+ Library 22 (21 repertoires)
+ Library 23 (35 repertoires)
+ Library 27 (38 repertoires)
+ Library 28 (41 repertoires)
+ Library 29 (38 repertoires)
+ Library 30 (39 repertoires)

## VDJPipe jobs

+ Library 13, 15, 16, 20, job id: 8c23efbf-d839-46e4-9c2e-24f3a150b4a7-007
+ Library 22, job id: 416577b3-3ea1-448c-9d67-8dd74a555f5f-007 (OLD, missing samples)
+ Library 22, job id: b2ca10f2-0fed-4a37-b9bb-719bbc2fe450-007 (with additional samples)
+ Library 23, job id: 1a7577ec-e35e-47d9-b24a-711d84ef7b39-007
+ Library 27, job id: 52bf3ece-4d58-4ae8-938b-78317613401a-007
+ Library 28, job id: d1d718cc-1622-429f-b30c-29721e0ecef1-007
+ Library 29, job id: 6ef52f66-e3c8-4113-b834-a456a15b25f3-007
+ Library 30, job id: 0429df15-0635-4791-a784-bb9b652cba44-007

For Library 30, it was discovered later that sample USCH028_S28_R1_001 was missing.
The vdjpipe job above is newer job with that sample included.

For Library 22, additional samples were added, so the whole library was re-processed.
The vdjpipe job above is newer job with that sample included.

+ Library 4, 6, 8, 9, 10, 23: these are some fairly old libraries so instead of
re-doing the pre-processing, I uploaded the post-filter FASTA file.

## IgBlast job

job id: 50524118-e81a-47f6-bc9d-2345969b9bfe-007
job id: adbaad15-6bbf-47fe-af6f-7553a44d9e21-007 (old, wrong germline DB)
job id: e1769e2b-1624-4d5f-83c7-0552e1b26877-007 (for missing USCH028_S28_R1_001 sample)


## RepCalc job

job id: d2296f11-4358-4972-9abd-263dde5c3287-007 (initial, not all repertoires, errors)

job id: 721704be-de75-44b6-b5de-1056e7f6a378-007 (not all repertoires, initial groups, errors)

These jobs have partial results. Did not analyze the errors for the first job as we needed
to add more repertoires and groups. For the second job, the error is sample
USCH028_S28_R1_001 is missing.

job id: a965a000-1436-499a-923f-a997816f7ac2-007 (not all repertoires, errors)

The IgBlast job (adbaad15-6bbf-47fe-af6f-7553a44d9e21-007) was run using the newer
AIRR germline DB, but unfortunately RepCalc doesn't support it yet so this job
has lots of errors. Also, there are some new samples to be added.
