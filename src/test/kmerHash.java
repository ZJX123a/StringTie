package test;

import java.io.IOException;
import test.create_graph;
import java.util.*;
import test.baseOptions;
import test.load_Read;

public class kmerHash {
	final int kmer_length = 25;
	final int min_exon_length = 20;
	final int fr_strand = 1;
	final int min_seed_cov=2;
	final float min_seed_entry=1.3f;
	final float g_min_ratio_non_error = 0.5f;
	final int pair_gap_length=200;
	final int max_pair_gap_length=500;
	final int min_anchor_length=5;
	final int min_function_count=2;
	final boolean is_paired_end=true;
	final float min_ratio_welds=0.04f;
	final int min_reads_span_junction=2;
	Map<Long, Vector<Integer>> kmer_hash = new HashMap<Long, Vector<Integer>>();
	
	Map<Long, Long> kmer_map = new HashMap<Long, Long>();
	List<Map.Entry<Long, Long>> list;

	public void readsToKmer(Map<Long, Vector<Integer>> kmer_hash) {
		// TODO Auto-generated method stub
		String temp_kmer;
		long intval = 0;
		for (int i = 0; i < load_Read.read_vector.size(); i++) {
			String read = load_Read.read_vector.get(i);
			for (int j = 0; j < read.length() - kmer_length + 1; j++) {
				temp_kmer = read.substring(j, j + kmer_length);
				intval = baseOptions.kmerToIntval(temp_kmer);
				Vector<Integer> readset = new Vector<Integer>();
				if (kmer_hash.containsKey(intval)) {
					readset = get_readset(intval, kmer_hash);
					readset.add(i);
				} else {
					readset.add(i);
				}
				kmer_hash.put(intval, readset);
			}

		}
	}

	public Vector<Integer> get_readset(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		return (Vector<Integer>) kmer_hash.get(intval);
	}

	private boolean ifexist(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	public long get_readset_count(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (ifexist(kmer_hash, intval)) {
			Vector readset = (Vector) kmer_hash.get(intval);
			long count = readset.size();
			return count;
		} else {
			return 0;
		}
	}

	private void delete_bad_kmers(Map<Long, Vector<Integer>> kmer_hash) {
		Vector delete_list = new Vector();
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
			Long key = iterator.next();
			Long candidate;
			candidates = get_forward_candidates(key, kmer_hash);
			long dominate_count = 0;
			if (candidates.size() > 0) {
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate = mapping.getKey();
					if (mapping.getValue() > 0) {
						long candidate_count = mapping.getValue();
						if (dominate_count == 0) {
							dominate_count = candidate_count;
						} else if ((float) candidate_count / dominate_count < g_min_ratio_non_error) {
							delete_list.add(candidate);
						}
					}
				}
			}
		}
		System.out.println("之前的hash");
		System.out.println(kmer_hash.size());
		if (!delete_list.isEmpty()) {
			for (int j = 0; j < delete_list.size(); j++) {
				long delete_kmer = (long) delete_list.get(j);
				kmer_hash.remove(delete_kmer);
			}
		}
		System.out.println(kmer_hash.size());
	}

	private boolean if_find_kmer(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	public List get_forward_candidates(long seed_kmer, Map<Long, Vector<Integer>> kmer_hash) {
		if (ifexist(kmer_hash, seed_kmer)) {
			Map forward_candidates = new HashMap();
			List<Map.Entry<Long, Long>> list_forward = new ArrayList<Map.Entry<Long, Long>>();
			long temp_intval = seed_kmer << (65 - kmer_length * 2) >> (63 - kmer_length * 2);
			long a[] = { 0l, 1l, 2l, 3l };
			for (int i = 0; i < 4; i++) {
				long temp2 = temp_intval;
				temp2 |= a[i];
				temp2 = baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
				if (ifexist(kmer_hash, temp2)) {
					long read_count = get_readset_count(kmer_hash, temp2);
					forward_candidates.put(temp2, read_count);
				}
				list_forward = new ArrayList<Map.Entry<Long, Long>>(forward_candidates.entrySet());

				// 通过比较器实现比较排序
				Collections.sort(list_forward, new Comparator<Map.Entry<Long, Long>>() {
					public int compare(Map.Entry<Long, Long> o1, Map.Entry<Long, Long> o2) {
						return o2.getValue().compareTo(o1.getValue()); // 倒序
					}
				});
			}
//			System.out.println(forward_candidates);
//			System.out.println(list_forward);
//			System.out.println();
			return list_forward;
		}
		return null;

	}

	public List get_reverse_candidates(long seed_kmer, Map<Long, Vector<Integer>> kmer_hash) {
		if (ifexist(kmer_hash, seed_kmer)) {
			Map reverse_candidates = new HashMap();
			List<Map.Entry<Long, Long>> list_reverse = new ArrayList<Map.Entry<Long, Long>>();
			long temp_intval = seed_kmer >> 2;
			long a[] = { 0l, 1l, 2l, 3l };
			for (int i = 0; i < 4; i++) {
				long temp2 = (a[i] << (2 * kmer_length - 2)) | temp_intval;
				temp2 = baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
				if (ifexist(kmer_hash, temp2)) {
					long read_count = get_readset_count(kmer_hash, temp2);
					//System.out.println(temp2+"    "+read_count);
					reverse_candidates.put(temp2, read_count);
				}

			}
			// if(reverse_candidates.size()==0){
			// return null;
			// }
			list_reverse = new ArrayList<Map.Entry<Long, Long>>(reverse_candidates.entrySet());

			// 通过比较器实现比较排序
			Collections.sort(list_reverse, new Comparator<Map.Entry<Long, Long>>() {
				public int compare(Map.Entry<Long, Long> o1, Map.Entry<Long, Long> o2) {
					return o2.getValue().compareTo(o1.getValue()); // 倒序
				}
			});
			return list_reverse;
		} else {
			return null;
		}

	}

	private void print_kmerhash(Map<Long, Vector<Integer>> kmer_hash) {
		System.out.println(kmer_hash.size());
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			Long key = iterator.next();
			System.out.println("key-----" + baseOptions.intvalToKmer(key, kmer_length) + "     " + key);
			System.out.println("value--------" + kmer_hash.get(key));
		}
	}

	public List sort_kmer(Map<Long, Vector<Integer>> kmer_hash) {

		Long key;
		Long count;
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			key = iterator.next();
			count = (long) kmer_hash.get(key).size();
			kmer_map.put(key, count);
		}
		list = new ArrayList<Map.Entry<Long, Long>>(kmer_map.entrySet());

		// 通过比较器实现比较排序
		Collections.sort(list, new Comparator<Map.Entry<Long, Long>>() {
			public int compare(Map.Entry<Long, Long> o1, Map.Entry<Long, Long> o2) {
				return o2.getValue().compareTo(o1.getValue()); // 倒序
			}
		});

		// for (Map.Entry<Long, Long> mapping : list) {
		// System.out.println(baseOptions.intvalToKmer(mapping.getKey(),
		// kmer_length) + ":" + mapping.getValue()
		// + " intval:" + mapping.getKey());
		// }
		return list;
	}

	public static void main(String args[]) throws IOException {
		load_Read.load_reads();
		kmerHash kh = new kmerHash();
		
		kh.readsToKmer(kh.kmer_hash);
		kh.sort_kmer(kh.kmer_hash);
		List listK = kh.sort_kmer(kh.kmer_hash);
		if (listK.size() == 0) {
			System.out.println("没有数据！");
			return;
		}
		else{
			System.out.println("有"+listK.size()+"个可用kmer");
		}
		int count=0;
		Set node_jihe=new HashSet();
		SplicingGraph sg = new SplicingGraph();
		System.out.println("初始used_kmers:"+sg.used_kmers.size());
		for(int i=0;i<listK.size();i++){
			if(!sg.has_been_used(kh.list.get(i).getKey())){
			    sg.init_trunk(kh,kh.list.get(i).getKey(),node_jihe,sg);
			    count++;
		    }
		}
		sg.rewrite_nodeSet(node_jihe);
		System.out.println("开始有"+sg.node_set.size()+"个点");
		sg.reunit_used_kmers(sg.node_set, kh);
		for(int i=0;i<sg.node_set.size();i++){
			sg.forward_check_and_extend(kh, i);
			sg.reverse_check_and_extend(kh, i);
		}
		sg.init_parents();
		for(int i=0;i<sg.node_set.size();i++){
			System.out.println("顶点编号："+i+"     顶点序列:"+sg.node_set.get(i).getSequence());
			System.out.println("父节点："+sg.node_set.get(i).getParents());
			System.out.println("子节点："+sg.node_set.get(i).getChildren());
		}
		System.out.println(sg.used_kmers.containsKey(49792356283881l));
	}
}
