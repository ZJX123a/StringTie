package test;
import java.io.IOException;
import java.util.*;
import test.baseOptions;
import test.load_Read;
public class kmerHash {
	static final int kmer_length=6;
	static Hashtable<Long, Vector<Integer>> kmer_hash=new Hashtable<Long, Vector<Integer>>();
	static float g_min_ratio_non_error = 0.5f;
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

	static void delete_bad_kmers() {
		int i = 0;
		Vector delete_list = new Vector();
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {

			Long key = iterator.next();
			System.out.println("处理第" + i + "个kmer" + "   " + baseOptions.intvalToKmer(key, kmer_length));
			Map candidates = get_forward_candidates(key);
			long dominate_count = 0;
			if (candidates.size() > 0) {
				for (Iterator<Long> its = candidates.keySet().iterator(); its.hasNext();) {
					Long candidate = its.next();
					if ((long) candidates.get(candidate) > 0) {
						long candidate_count = (long) candidates.get(candidate);
						System.out.println("其替补kmer是：" + baseOptions.intvalToKmer(candidate, kmer_length) + "   count:"
								+ candidate_count);
						if (dominate_count == 0) {
							dominate_count = candidate_count;
						} else if ((float) candidate_count / dominate_count < g_min_ratio_non_error) {
							delete_list.add(candidate);
						}
					}
				}
			}
			i++;
			System.out.println();
		}
		System.out.println("之前的hash");
		print_kmerhash();
		System.out.println();
		if (!delete_list.isEmpty()) {
			for (int j = 0; j < delete_list.size(); j++) {
				long delete_kmer = (long) delete_list.get(j);
				kmer_hash.remove(delete_kmer);
			}
		}
		System.out.println("之后的hash");
		System.out.println(kmer_hash.size());
		print_kmerhash();

	}
   static boolean if_find_kmer(long intval){
	   if(kmer_hash.containsKey(intval)){
		   return true;
	   }
	   else{
		   return false;
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
		    	temp2=baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
		    	if(ifexist(kmer_hash,temp2)){
			      long read_count=get_readset_count(kmer_hash,temp2);
			      forward_candidates.put(temp2, read_count);
		    	}
			   // System.out.println(temp2);
			    String str=baseOptions.intvalToKmer(temp2,kmer_length);
			    //System.out.println(str);
		    }
		  //  System.out.println(forward_candidates);
		    return forward_candidates;
	   }
	return null;
	   
   }
   static void print_kmerhash(){
	   System.out.println(kmer_hash.size());
	   for(Iterator<Long> iterator=kmer_hash.keySet().iterator();iterator.hasNext();){
		   Long key=iterator.next();
		   System.out.println("key-----"+baseOptions.intvalToKmer(key,kmer_length));
		 //  System.out.println("value--------"+kmer_hash.get(key));
		   }
   }
   public static void main(String args[]) throws IOException{
	   load_Read.load_reads();
	   readsToKmer();
	//   print_kmerhash();
	   delete_bad_kmers();
	  // long seed_int=baseOptions.kmerToIntval("AGAAACCAACAGTGTGCTTTTAATA");
	   //System.out.println(seed_int);
	 //  get_forward_candidates(baseOptions.kmerToIntval("TCCCTT"));
	 //  delete_bad_kmers();
   }
}
