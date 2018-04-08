package test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import test.kmerHash;

public class create_graph {
	Map<Long, Long> used_kmers = new HashMap<Long, Long>();
	private  String forward_extend(long kmer_int,Vector<Long> bifurcation,kmerHash kh){
		List<Map.Entry<Long,Long>> candidates=new ArrayList<Map.Entry<Long,Long>>();
		boolean flag = false;
		long sum=0l;
		String kmer_string=baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while(true){
			candidates=kh.get_forward_candidates(kmer_int,kh.kmer_hash);
			//System.out.println("����kmer:"+baseOptions.intvalToKmer(kmer_int, kh.kmer_length));
			if(candidates.size()==0){
				break;
			}
			int count=0;
			long candidate = 0l;
			long cov=0l;
			//System.out.println("��candidatesΪ��");
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate=mapping.getKey();
				//System.out.println("candidate:"+baseOptions.intvalToKmer(candidate, kh.kmer_length)+"   read_count:"+mapping.getValue());
				if(!has_been_used(candidate)){//�����δ��ʹ�õ�candidate
				//	System.out.println("ѡ�У�"+baseOptions.intvalToKmer(candidate, kh.kmer_length));
					flag=true;
					cov=mapping.getValue();
					//�����������δ��ʹ�õ�candidate
					if(count<candidates.size()-1){
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
			//	System.out.println();
			}
			if(flag==true){
				used_kmers.put(candidate, cov);
				sum+=cov;
				int base_last_int=(int) (candidate&3l);
				char base_last_char=baseOptions.intToBase(base_last_int);
				kmer_string+=base_last_char;
				kmer_int=candidate;
			}
		}
		return kmer_string;
	}

	//������չ
	private  String reverse_extend(long kmer_int,Vector<Long> bifurcation,kmerHash kh){
		List<Map.Entry<Long,Long>> candidates=new ArrayList<Map.Entry<Long,Long>>();
		boolean flag = false;
		long sum=0l;
		String kmer_string=baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while(true){
			candidates=kh.get_reverse_candidates(kmer_int,kh.kmer_hash);
			System.out.println("����kmer:"+baseOptions.intvalToKmer(kmer_int, kh.kmer_length));
			if(candidates.size()==0){
				break;
			}
			int count=0;
			long candidate = 0l;
			long cov=0l;
			System.out.println("��candidatesΪ��");
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate=mapping.getKey();
				System.out.println("candidate:"+baseOptions.intvalToKmer(candidate, kh.kmer_length)+"   read_count:"+mapping.getValue());
				if(!has_been_used(candidate)){//�����δ��ʹ�õ�candidate
					System.out.println("ѡ�У�"+baseOptions.intvalToKmer(candidate, kh.kmer_length));
					flag=true;
					cov=mapping.getValue();
					//�����������δ��ʹ�õ�candidate
					if(count<candidates.size()-1){
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
				System.out.println();
			}
			if(flag==true){
				used_kmers.put(candidate, cov);
				sum+=cov;
				int base_first_int=(int) (candidate>>(2*kh.kmer_length-2));
				char base_first_char=baseOptions.intToBase(base_first_int);
				kmer_string=base_first_char+kmer_string;
				kmer_int=candidate;
			}
		}
		return kmer_string;
	}

	
	private  boolean has_been_used(Long candidate) {
		if (used_kmers.containsKey(candidate)) {
			return true;
		} else {
			return false;
		}
	}

	public void init_graph(kmerHash kh) {
		// TODO Auto-generated method stub
//		System.out.println(kh.list.get(0));
//		System.out.println("��ʼkmer:"+baseOptions.intvalToKmer(kh.list.get(0).getKey(), kh.kmer_length));
//		System.out.println("**************");
//		used_kmers.put(kh.list.get(0).getKey(), kh.list.get(0).getValue());
//		System.out.println(used_kmers.toString());
//		System.out.println("������");
		Vector<Long> bifurcation=new Vector<Long>();
		String right=forward_extend(kh.list.get(0).getKey(),bifurcation,kh);
		String left=reverse_extend(kh.list.get(0).getKey(),bifurcation,kh);
		String trunk=left+right.substring(kh.kmer_length);
	//	System.out.println("�����չ��"+left);
	//	System.out.println("�ұ���չ��"+right);
		System.out.println("******************"+trunk);
		System.out.println();
	}

}
