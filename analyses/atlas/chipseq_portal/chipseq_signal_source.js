import {BigWig} from '@gmod/bbi';
import {BlobFile} from 'generic-filehandle2';
export async function openBigWig(blob,chroms){
  const reader=new BigWig({filehandle:new BlobFile(blob)}),header=await reader.getHeader();
  const refs=Object.values(header.refsByNumber);
  if(!refs.length||refs.some(ref=>!Object.hasOwn(chroms,ref.name)||chroms[ref.name]!==ref.length))throw Error('This signal file does not match the SM_V10 chromosome names and lengths.');
  return {features:async(locus,pixels)=>reader.getFeatures(locus.chrom,locus.start,locus.end,{basesPerSpan:(locus.end-locus.start)/pixels})};
}
