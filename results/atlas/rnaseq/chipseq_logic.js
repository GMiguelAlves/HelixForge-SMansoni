(function(root,factory){const api=factory();if(typeof module==='object')module.exports=api;else root.ChipLogic=api;})(globalThis,()=>{
  const overlap=(a,b)=>a.chrom===b.chrom&&a.start<b.end&&a.end>b.start;
  function chromosome(input,chroms){
    const text=String(input).trim().toUpperCase();
    if(Object.hasOwn(chroms,text))return text;
    const alias=text.replace(/^CHR/,'').replace(/^SM_V10_/,'');
    const name='SM_V10_'+({M:'MITO',MT:'MITO',W:'WSR'}[alias]||alias);
    return Object.hasOwn(chroms,name)?name:null;
  }
  function parseLocus(input,chroms){
    const match=String(input).trim().match(/^([^:\s]+):\s*([\d,]+)\s*-\s*([\d,]+)$/);
    if(!match)return null;
    const chrom=chromosome(match[1],chroms),start=Number(match[2].replaceAll(',','')),end=Number(match[3].replaceAll(',',''));
    if(!chrom||!Number.isSafeInteger(start)||!Number.isSafeInteger(end)||start<1||end<start||end>chroms[chrom])throw Error('Use a valid SM_V10 interval: chromosome:start-end (1-based, inclusive).');
    return {chrom,start:start-1,end};
  }
  function clamp(locus,chroms){const max=chroms[locus.chrom];const width=Math.max(20,Math.min(max,Math.round(locus.end-locus.start)));const start=Math.max(0,Math.min(max-width,Math.round(locus.start)));return {chrom:locus.chrom,start,end:start+width};}
  function move(locus,fraction,chroms){const d=(locus.end-locus.start)*fraction;return clamp({...locus,start:locus.start+d,end:locus.end+d},chroms);}
  function zoom(locus,factor,chroms){const mid=(locus.start+locus.end)/2,width=(locus.end-locus.start)*factor;return clamp({...locus,start:mid-width/2,end:mid+width/2},chroms);}
  function geneMatches(genes,query){const q=query.toLowerCase().trim();if(!q)return [];return genes.filter(g=>[g.id,g.description,g.previous].some(v=>String(v).toLowerCase().includes(q))).sort((a,b)=>Number(b.id.toLowerCase()===q)-Number(a.id.toLowerCase()===q)||a.id.localeCompare(b.id));}
  function matches(row,filters){return (!filters.study||row.study===filters.study)&&(!filters.stage||row.stage===filters.stage)&&(!filters.factor||row.factor===filters.factor)&&(!filters.category||row.category===filters.category)&&(!filters.chrom||row.chrom===filters.chrom)&&(!filters.locus||overlap(row,filters.locus))&&(row.kind==='library'||(!filters.minSupport||row.support>=Number(filters.minSupport))&&(!filters.ambiguity||row.ambiguity===filters.ambiguity));}
  function pack(rows,limit=10){const ends=[],placed=[];for(const row of [...rows].sort((a,b)=>a.start-b.start||a.end-b.end)){let lane=ends.findIndex(end=>end<=row.start);if(lane<0)lane=ends.length;if(lane>=limit)continue;ends[lane]=row.end;placed.push({...row,lane});}return {rows:placed,lanes:Math.max(1,ends.length),omitted:rows.length-placed.length};}
  function bins(features,locus,count=900){const out=Array(count).fill(null),weights=Array(count).fill(0),span=locus.end-locus.start;for(const f of features){if(!Number.isFinite(f.score)||f.end<=locus.start||f.start>=locus.end)continue;const a=Math.max(0,Math.floor((f.start-locus.start)/span*count)),b=Math.min(count-1,Math.floor((f.end-1-locus.start)/span*count));for(let i=a;i<=b;i++){const width=Math.max(0,Math.min(f.end,locus.start+(i+1)*span/count)-Math.max(f.start,locus.start+i*span/count));if(width>0){out[i]=(out[i]??0)+f.score*width;weights[i]+=width;}}}return out.map((v,i)=>weights[i]?v/weights[i]:null);}
  const locusLabel=l=>`${l.chrom}:${(l.start+1).toLocaleString('en-US')}-${l.end.toLocaleString('en-US')}`;
  return {overlap,chromosome,parseLocus,clamp,move,zoom,geneMatches,matches,pack,bins,locusLabel};
});
