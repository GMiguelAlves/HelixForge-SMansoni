// Upstream Volvox parser fixture: never included as S. mansoni atlas evidence.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import {openBigWig} from '../analyses/atlas/chipseq_portal/vendor/bbi_reader.js';
const fixture=new Blob([fs.readFileSync(new URL('./fixtures/chipseq_parser_volvox.bw',import.meta.url))]);
const reader=await openBigWig(fixture,{ctgA:50001});
const features=await reader.features({chrom:'ctgA',start:0,end:50001},800);
assert(features.length>0);
assert(features.every(f=>f.end>f.start&&Number.isFinite(f.score)));
await assert.rejects(()=>openBigWig(fixture,{SM_V10_1:87984048}),/does not match/);
await assert.rejects(()=>openBigWig(fixture,{ctgA:50000}),/does not match/);
await assert.rejects(()=>openBigWig(new Blob(['invalid BigWig']),{ctgA:50001}));
console.log(`PASS: browser BigWig reader decoded ${features.length} features; incompatible chromosome names, lengths and malformed files rejected.`);
