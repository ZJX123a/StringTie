package test;
import java.io.IOException;
import java.util.*;
import test.baseOptions;
import test.load_Read;
public class kmerHash {
	static final int kmer_length=25;
	static Hashtable<Long, Vector<Integer>> kmer_hash=new Hashtable<Long, Vector<Integer>>();
	public static void readsToKmer() {
		// TODO Auto-generated method stub
		String temp_kmer;
		long intval = 0;
		for( int i=0;i<load_Read.read_vector.size();i++){
			String read=load_Read.read_vector.get(i);
			for(int j=0;j<read.length()-kmer_length+1;j++){
				temp_kmer=read.substring(j,j+kmer_length);
				intval=baseOptions.kmerToIntval(temp_kmer);
				//System.out.println("j::"+j+"     kl::"+kmer_length+"       "+temp_kmer+"     "+intval);
				Vector<Integer> readset=new Vector<Integer>();
				//readset.add(i);
				if(kmer_hash.containsKey(intval)){
					readset=get_readset(intval,kmer_hash);
					readset.add(i);
				}
				else{
					readset.add(i);
				}
				kmer_hash.put(intval, readset);
			}
			//System.out.println(kmer_hash.get(intval));
			
		}
		System.out.println(kmer_hash.toString());
	}
	
	static Vector<Integer> get_readset(long intval,Hashtable<Long, Vector<Integer>> kmer_hash){
		return (Vector<Integer>) kmer_hash.get(intval);
	}
	static boolean ifexist(Hashtable kmer_hash,long intval){
		if(kmer_hash.containsKey(intval)){
			return true;
		}
		else{
			return false;
		}
	}
   static long get_readset_count(Hashtable kmer_hash,long intval){
    	if(ifexist(kmer_hash,intval)){
    		Vector readset=(Vector) kmer_hash.get(intval);
    		long count=readset.size();
    		return count;
    	}
    	else{
    		return 0;
    	}
    }
   static void delete_bad_kmers(){
	   for(Iterator<Long> iterator=kmer_hash.keySet().iterator();iterator.hasNext();){
		   Long key=iterator.next();
		   System.out.println("key-----"+key);
		   System.out.println("value--------"+kmer_hash.get(key));
		   }
   }
   static Map get_forward_candidates(long seed_kmer){
	   if(ifexist(kmer_hash,seed_kmer)){
		    Map forward_candidates=new HashMap();
		    long temp_intval=seed_kmer<<(65-kmer_length*2)>>(63-kmer_length*2);
		    int a[]={0,1,2,3};
		    for(int i=0;i<4;i++){
		    	long temp2=temp_intval;
		    	temp2|=a[i];
			    long read_count=get_readset_count(kmer_hash,temp2);
			    forward_candidates.put(temp2, read_count);
			    System.out.println(temp2);
			    String str=baseOptions.intvalToKmer(temp2,kmer_length);
			    System.out.println(str);
		    }
		    
		    return forward_candidates;
	   }
	return null;
	   
   }
   public static void main(String args[]) throws IOException{
	   load_Read.load_reads();
	   readsToKmer();
	  // delete_bad_kmers();
	   long seed_int=baseOptions.kmerToIntval("AGAAACCAACAGTGTGCTTTTAATA");
	   System.out.println(seed_int);
	   get_forward_candidates(seed_int);
   }
}
