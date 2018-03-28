package test;
import java.util.*;
import test.baseOptions;
public class kmerHash {
	static final int kmer_length=6;
	Hashtable kmer_hash=new Hashtable();
	public void readsToKmer() {
		// TODO Auto-generated method stub
		String temp_kmer;
		for( int i=0;i<load_Read.read_vector.size();i++){
			String read=load_Read.read_vector.get(i);
			for(int j=0;j<read.length()-kmer_length+1;j++){
				temp_kmer=read.substring(j,j+kmer_length);
				long intval=baseOptions.kmerToIntval(temp_kmer);
				System.out.println("j::"+j+"     kl::"+kmer_length+"       "+temp_kmer+"     "+intval);
				Vector readset=new Vector();
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
			System.out.println("第"+i+"条数据：：：：：");
			System.out.println(kmer_hash.toString());
		}
	}
	
	Vector get_readset(long intval,Hashtable kmer_hash){
		
		return (Vector) kmer_hash.get(intval);
	}

}
