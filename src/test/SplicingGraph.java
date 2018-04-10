package test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import test.SplicingGraph.Node;

public class SplicingGraph {
	Map<Long, Long> used_kmers = new HashMap<Long, Long>();
	Vector<Node> node_set = new Vector<Node>();
	Set<Long> restore_kmers = new HashSet<Long>();
	Set<Integer> forward_branches= new HashSet<Integer>();
	public class Node {
		private Vector<Integer> parents = new Vector<Integer>();
		private Vector<Integer> children = new Vector<Integer>();
		private String sequence;

		public boolean addParents(int parent) {
			if (parent < 0) {
				return false;
			} else {
				if (parents.size() != 0) {
					for (int i = 0; i < parents.size(); ++i) {
						if (parents.get(i) == parent) // if exist already
							return false;
					}
				}
				this.parents.add(parent);
				return true;
			}
		}

		public List getParents() {
			return parents;
		}

		public boolean addChildren(int child) {
			if (child < 0) {
				return false;
			} else {
				if (children.size() != 0) {
					for (int i = 0; i < children.size(); ++i) {
						if (children.get(i) == child) // if exist already
							return false;
					}
				}
				this.children.add(child);
				return true;
			}
		}

		public Vector getChildren() {
			return children;
		}

		public void setSequence(String seq) {
			sequence = seq;
		}

		public String getSequence() {
			return sequence;
		}

		public boolean isChild(int child) {
			if (child < 0) {
				return false;
			} else {
				for (int i = 0; i < children.size(); i++) {
					if (children.get(i) == child) {
						return true;
					}
				}
				return false;
			}
		}

		public boolean isParent(int parent) {
			if (parent < 0) {
				return false;
			} else {
				for (int i = 0; i < parents.size(); i++) {
					if (parents.get(i) == parent) {
						return true;
					}
				}
				return false;
			}
		}

		public boolean deleteChild(long child) {
			// TODO Auto-generated method stub
			if (child < 0) {
				return false;
			} else {
				for (int i = 0; i < children.size(); i++) {
					if (children.get(i) == child) {
						children.remove(i);
						return true;
					}
				}
				return false;
			}
		}

		public boolean deleteParent(long parent) {
			// TODO Auto-generated method stub
			if (parent < 0) {
				return false;
			} else {
				for (int i = 0; i < parents.size(); i++) {
					if (parents.get(i) == parent) {
						parents.remove(i);
						return true;
					}
				}
				return false;
			}
		}

		public void clearChildren() {
			// TODO Auto-generated method stub
			children.clear();
		}

		public void clearParents() {
			parents.clear();
		}

		public void clearAll() {
			// TODO Auto-generated method stub
			children.clear();
			parents.clear();
			sequence = "";
		}
	}

	private String forward_extend(long kmer_int, Vector<Long> bifurcation, kmerHash kh) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		boolean flag = false;
		long sum = 0l;
		String kmer_string = baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while (true) {
			candidates = kh.get_forward_candidates(kmer_int, kh.kmer_hash);
			// System.out.println("处理kmer:"+baseOptions.intvalToKmer(kmer_int,
			// kh.kmer_length));
			if (candidates.size() == 0) {
				break;
			}
			int count = 0;
			long candidate = 0l;
			long cov = 0l;
			// System.out.println("其candidates为：");
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				// System.out.println("candidate:"+baseOptions.intvalToKmer(candidate,
				// kh.kmer_length)+" read_count:"+mapping.getValue());
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					// System.out.println("选中："+baseOptions.intvalToKmer(candidate,
					// kh.kmer_length));
					flag = true;
					cov = mapping.getValue();
					// 如果还有其他未被使用的candidate
					if (count < candidates.size() - 1) {
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
				// System.out.println();
			}
			if (flag == true) {
				used_kmers.put(candidate, cov);
				sum += cov;
				int base_last_int = (int) (candidate & 3l);
				char base_last_char = baseOptions.intToBase(base_last_int);
				kmer_string += base_last_char;
				kmer_int = candidate;
			}
		}
		return kmer_string;
	}

	// 反向扩展
	private String reverse_extend(long kmer_int, Vector<Long> bifurcation, kmerHash kh) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		boolean flag = false;
		long sum = 0l;
		String kmer_string = baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while (true) {
			candidates = kh.get_reverse_candidates(kmer_int, kh.kmer_hash);
			// System.out.println("处理kmer:" + baseOptions.intvalToKmer(kmer_int,
			// kh.kmer_length));
			if (candidates.size() == 0) {
				break;
			}
			int count = 0;
			long candidate = 0l;
			long cov = 0l;
			// System.out.println("其candidates为：");
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				// System.out.println("candidate:" +
				// baseOptions.intvalToKmer(candidate, kh.kmer_length) + "
				// read_count:"
				// + mapping.getValue());
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					// System.out.println("选中：" +
					// baseOptions.intvalToKmer(candidate, kh.kmer_length));
					flag = true;
					cov = mapping.getValue();
					// 如果还有其他未被使用的candidate
					if (count < candidates.size() - 1) {
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
				System.out.println();
			}
			if (flag == true) {
				used_kmers.put(candidate, cov);
				sum += cov;
				int base_first_int = (int) (candidate >> (2 * kh.kmer_length - 2));
				char base_first_char = baseOptions.intToBase(base_first_int);
				kmer_string = base_first_char + kmer_string;
				kmer_int = candidate;
			}
		}
		return kmer_string;
	}

	public boolean has_been_used(Long candidate) {
		if (used_kmers.containsKey(candidate)) {
			return true;
		} else {
			return false;
		}
	}

	public void extend_again(kmerHash kh, Vector<Long> bifurcation) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		if (node_set.get(0).sequence.length() > 2 * kh.kmer_length + kh.kmer_length) {
			for (int i = 0; i < 2 * kh.kmer_length; i++) {
				String kmer = node_set.get(0).sequence.substring(i, kh.kmer_length + i);
				// System.out.println(kmer + " &");
				long intval = baseOptions.kmerToIntval(kmer);
				candidates = kh.get_reverse_candidates(intval, kh.kmer_hash);
				// System.out.println("candidates.size():" + candidates.size());
				if (candidates.size() <= 1) {
					continue;
				}
				boolean flag = false;
				long candidate;
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate = mapping.getKey();
					// System.out.println("candidate:" + candidate);
					if (!has_been_used(candidate)) {
						flag = true;
						break;
					}
				}
				if (flag == true) {
					String sequence = reverse_extend(intval, bifurcation, kh);
					if (sequence.length() > kh.kmer_length + kh.min_exon_length) {
						node_set.get(0).sequence = sequence + node_set.get(0).sequence.substring(i + kh.kmer_length);
						break;
					}
				}
			}
		}
		if (node_set.get(0).sequence.length() > 3 * kh.kmer_length + kh.kmer_length) {
			for (int i = node_set.get(0).sequence.length() - 3 * kh.kmer_length; i < node_set.get(0).sequence.length()
					- kh.kmer_length; i++) {
				String kmer = node_set.get(0).sequence.substring(i, kh.kmer_length + i);
				long intval = baseOptions.kmerToIntval(kmer);
				candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
				// System.out.println("forward_candidates:"+candidates.size());
				if (candidates == null) {
					continue;
				}
				boolean flag = false;
				long candidate;
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate = mapping.getKey();
					if (!has_been_used(candidate)) {
						flag = true;
						break;
					}
				}
				if (flag == true) {
					String sequence = forward_extend(intval, bifurcation, kh);
					if (sequence.length() > kh.kmer_length + kh.min_exon_length) {
						node_set.get(0).sequence = node_set.get(0).sequence.substring(0, i) + sequence;
						break;
					}
				}
			}
		}
	}

	public int add_node(Node node) {
		node_set.add(node);
		System.out.println("节点个数：" + node_set.get(0));
		return node_set.size() - 1;

	}

	public void forward_extend_use_pairInfo(kmerHash kh, Vector<String> data, int node_index,
			Vector<Long> bifurcation) {
		// 利用paired_end_reads进行扩展 node_index为顶点标号
		int middle_read_id = data.size() / 2;
		long extend_val = 0l;
		while (true) {
			int start_pos;
			int length = node_set.get(node_index).sequence.length();
			if (length > 2 * kh.kmer_length) {
				start_pos = length - 2 * kh.kmer_length;
			} else {
				start_pos = 0;
			}
			String check = node_set.get(node_index).sequence.substring(start_pos);
			boolean extend_flag = false;
			Set<Integer> reads = new HashSet<Integer>();
			set_reads(kh, check, reads);
			if (reads.size() == 0) {
				return;
			}
			Iterator<Integer> it = reads.iterator();
			int paired_end_read_id = 0;
			String extend_str = "";
			while (it.hasNext()) {
				int read_id = it.next();
				if (kh.fr_strand == 1) { // 2-> <-1 此时应该把1反过来 根据2找1
					if (read_id < middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id - middle_read_id;
					}
				} else if (kh.fr_strand == 2) {
					if (read_id >= middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id + middle_read_id;
					}
				}
				String paired_end_read = load_Read.read_vector.get(paired_end_read_id);
				long max = 0l;
				long max_read_kmer = 0l;
				for (int i = 0; i < paired_end_read.length() - kh.kmer_length; i++) {
					String kmer = paired_end_read.substring(i, i + kh.kmer_length);
					long intval = baseOptions.kmerToIntval(kmer);
					if (has_been_used(intval)) {
						max = 0;
						break;
					} else if (kh.get_readset_count(kh.kmer_hash, intval) > max) {
						max = kh.get_readset_count(kh.kmer_hash, intval);
						max_read_kmer = intval;
					}
				}
				// System.out.println("entry:"+baseOptions
				// .computeEntropy(baseOptions.intvalToKmer(max_read_kmer,
				// kh.kmer_length)));
				if (max >= kh.min_seed_cov && baseOptions
						.computeEntropy(baseOptions.intvalToKmer(max_read_kmer, kh.kmer_length)) > kh.min_seed_entry) {
					long stop_kmer = baseOptions
							.kmerToIntval(node_set.get(node_index).sequence.substring(length - kh.kmer_length));
					// String str = "";
					StringBuffer str = new StringBuffer("");
					extend_val = find_head_kmer(kh, max_read_kmer, stop_kmer, str);
					// System.out.println("str:" + str);
					String extend_kmer = baseOptions.intvalToKmer(extend_val, kh.kmer_length);
					String anchor = extend_kmer.substring(0, 5);
					int start = check.indexOf(anchor);
					if (start != -1) { // 如果找到了
						if (is_similar(check.substring(start), extend_kmer, 'F')) {
							node_set.get(node_index).sequence = node_set.get(node_index).sequence.substring(0,
									start_pos + start) + str;
							extend_flag = true;
						}
					} else {
						if (((int) str.length() > kh.pair_gap_length - 80) && (str.length() > extend_str.length()))
							extend_str = str.toString();
					}
					if (extend_flag) {
						add_used_kmers(kh, str.toString());
						break;
					}
				}

			}

			if (extend_flag) { // 那么将得到的str再向前扩展
				String kmer = node_set.get(node_index).sequence
						.substring(node_set.get(node_index).sequence.length() - kh.kmer_length);
				long kmer_intval = baseOptions.kmerToIntval(kmer);
				String str = forward_extend(kmer_intval, bifurcation, kh);
				node_set.get(node_index).sequence = node_set.get(node_index).sequence + str.substring(kh.kmer_length);
				if (str.length() < 2 * kh.kmer_length)
					return;
			} else {

				if (((int) extend_str.length() > kh.pair_gap_length - 80)
						&& ((int) extend_str.length() < kh.max_pair_gap_length)) {
					node_set.get(node_index).sequence = node_set.get(node_index).sequence + extend_str;
					add_used_kmers(kh, extend_str);
				}

				return;
			}
		}
	}

	public void reverse_extend_use_pairInfo(kmerHash kh, Vector<String> data, int node_index,
			Vector<Long> bifurcation) {
		int middle_read_id = data.size() / 2;
		while (true) {
			String check = node_set.get(node_index).sequence.substring(0, 2 * kh.kmer_length);
			Set<Integer> reads = new HashSet<Integer>();
			set_reads(kh, check, reads);
			if (reads.size() == 0) {
				return;
			}
			long extend_val = 0l;
			boolean extend_flag = false;
			Iterator<Integer> it = reads.iterator();
			int paired_end_read_id = 0;
			String extend_str = "";
			while (it.hasNext()) {
				int read_id = it.next();
				if (kh.fr_strand == 1) { // --2--> ..........
											// <--1--//反向扩展，应该根据1找2，连接2
					if (read_id >= middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id + middle_read_id;
					}
				} else if (kh.fr_strand == 2) { // --1--> .......... <--2--
					if (read_id < middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id - middle_read_id;
					}
				}
				String paired_end_read = load_Read.read_vector.get(paired_end_read_id);
				long max = 0l;
				long max_read_kmer = 0l;
				for (int i = 0; i < paired_end_read.length() - kh.kmer_length; i++) {
					String kmer = paired_end_read.substring(i, i + kh.kmer_length);
					long intval = baseOptions.kmerToIntval(kmer);
					if (has_been_used(intval)) {
						max = 0;
						break;
					} else if (kh.get_readset_count(kh.kmer_hash, intval) > max) {
						max = kh.get_readset_count(kh.kmer_hash, intval);
						max_read_kmer = intval;
					}
				}

				if (max >= kh.min_seed_cov && baseOptions
						.computeEntropy(baseOptions.intvalToKmer(max_read_kmer, kh.kmer_length)) > kh.min_seed_entry) {
					long stop_kmer = baseOptions
							.kmerToIntval(node_set.get(node_index).sequence.substring(0, kh.kmer_length));
					StringBuffer str = new StringBuffer("");
					extend_val = find_tail_kmer(kh, max_read_kmer, stop_kmer, str);
					// System.out.println("str:" + str);
					String extend_kmer = baseOptions.intvalToKmer(extend_val, kh.kmer_length);
					String anchor = extend_kmer.substring(kh.kmer_length - 5);
					int start = check.indexOf(anchor);
					if (start != -1) { // 如果找到了
						if (is_similar(check.substring(0, start + 5), extend_kmer, 'R')) {
							node_set.get(node_index).sequence = str.toString()
									+ node_set.get(node_index).sequence.substring(start + 5);
							extend_flag = true;
						}
					} else {
						if (((int) str.length() > kh.pair_gap_length - 80) && (str.length() > extend_str.length()))
							extend_str = str.toString();
					}
					if (extend_flag) {
						add_used_kmers(kh, str.toString());
						break;
					}
				}
			}
			if (extend_flag) { // 那么将得到的str再逆向扩展
				String kmer = node_set.get(node_index).sequence.substring(0, kh.kmer_length);
				long kmer_intval = baseOptions.kmerToIntval(kmer);
				String str = reverse_extend(kmer_intval, bifurcation, kh);
				node_set.get(node_index).sequence = str.substring(0, str.length() - kh.kmer_length)
						+ node_set.get(node_index).sequence;
				if (str.length() < 2 * kh.kmer_length)
					return;
			} else {

				if (((int) extend_str.length() > kh.pair_gap_length - 80)
						&& ((int) extend_str.length() < kh.max_pair_gap_length)) {
					node_set.get(node_index).sequence = extend_str + node_set.get(node_index).sequence;
					add_used_kmers(kh, extend_str);
				}

				return;
			}

		}
	}

	public long find_tail_kmer(kmerHash kh, long seed_intval, long stop_kmer, StringBuffer str) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		long intval = seed_intval;
		str.append(baseOptions.intvalToKmer(intval, kh.kmer_length));
		Map<Long, Boolean> use_kmers = new HashMap<Long, Boolean>();
		while (true) {
			candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
			if (candidates.size() == 0) {
				break;
			}
			long candidate = 0l;
			boolean flag = false;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!use_kmers.containsKey(candidate)) {
					flag = true;
					break;
				}
			}
			// 如果所有的candidate都在use_kmers中，或者candidate为stop_kmer 退出while循环
			if ((!flag) || candidate == stop_kmer) {
				break;
			}
			use_kmers.put(candidate, true);
			intval = candidate;
			int base_int = (int) (candidate & 3l);
			char base_char = baseOptions.intToBase(base_int);
			str.append(base_char);
		}
		// System.out.println("real_str:" + str);
		return intval;
	}

	public boolean compatible(String seq1, String seq2, kmerHash kh) {
		for (int i = 0; i < seq2.length() - kh.kmer_length; i++) {
			String kmer = seq2.substring(i, i + kh.kmer_length);
			int index = seq1.indexOf(kmer);
			if (index != -1) {
				if (index > i) {
					return is_aligned(seq1.substring(index - i), seq2);
				} else {
					return is_aligned(seq1, seq2.substring(i - index));
				}
			}

		}
		return false;

	}

	public boolean is_aligned(String seq1, String seq2) {
		int mismatch = 0;
		int length;
		if (seq1.length() >= seq2.length()) {
			length = seq2.length();
		} else {
			length = seq1.length();
		}
		for (int i = 0; i < length; ++i) {
			if (seq1.charAt(i) != seq2.charAt(i))
				mismatch++;
		}
		return (mismatch <= 2);
	}

	public boolean add_used_kmers(kmerHash kh, String str) {
		int length = str.length();
		if (length >= kh.kmer_length) {
			for (int i = 0; i <= length - kh.kmer_length; ++i) {
				String kmer = str.substring(i, i + kh.kmer_length);
				long intval = baseOptions.kmerToIntval(kmer);
				long cov = kh.get_readset_count(kh.kmer_hash, intval);
				used_kmers.put(intval, cov);
			}
		}

		return true;
	}

	public boolean is_similar(String seq1, String seq2, char direction) {
		int mismatch = 0;
		int length;
		if (seq1.length() >= seq2.length()) {
			length = seq2.length();
		} else {
			length = seq1.length();
		}
		if (length == 0) {
			if (seq1.length() == 0 && seq2.length() == 0) {
				return true;
			} else {
				return false;
			}
		}
		if (direction == 'F') { // 从前往后检查
			for (int i = 0; i < length; ++i) {
				if (seq1.charAt(i) != seq2.charAt(i))
					mismatch++;
			}
		} else {
			for (int i = 0; i < length; ++i) {
				if (seq1.charAt(seq1.length() - i - 1) != seq2.charAt(seq2.length() - i - 1))
					mismatch++;
			}
		}
		if ((float) mismatch / length < 0.35) {
			return true;
		} else {
			return false;
		}
	}

	public long find_head_kmer(kmerHash kh, long seed_intval, long stop_kmer, StringBuffer str) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		long intval = seed_intval;
		str.append(baseOptions.intvalToKmer(intval, kh.kmer_length));
		Map<Long, Boolean> use_kmers = new HashMap<Long, Boolean>();
		while (true) {
			candidates = kh.get_reverse_candidates(intval, kh.kmer_hash);
			if (candidates.size() == 0) {
				break;
			}
			long candidate = 0l;
			boolean flag = false;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!use_kmers.containsKey(candidate)) {
					flag = true;
					break;
				}
			}
			// 如果所有的candidate都在use_kmers中，或者candidate为stop_kmer 退出while循环
			if ((!flag) || candidate == stop_kmer) {
				break;
			}
			use_kmers.put(candidate, true);
			intval = candidate;
			int base_int = (int) ((candidate >> (kh.kmer_length * 2 - 2)) & 3l);
			char base_char = baseOptions.intToBase(base_int);
			str.insert(0, base_char);
		}
		// System.out.println("real_str:" + str);
		return intval;

	}

	public void set_reads(kmerHash kh, String seq, Set reads) {
		// 返回seq覆盖到的所有reads的编号 不重复
		if (seq.length() < kh.kmer_length) {
			return;
		}
		for (int i = 0; i < seq.length() - kh.kmer_length; i++) {
			String kmer = seq.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			Vector read_set = kh.get_readset(intval, kh.kmer_hash);
			if (read_set.size() == 0) {
				continue;
			}
			for (int j = 0; j < read_set.size(); j++) {
				reads.add(read_set.get(j));
			}
		}
	}

	public boolean init_trunk(kmerHash kh, long seed_val, Set node_jihe) {
		// TODO Auto-generated method stub
		used_kmers.put(seed_val, kh.get_readset_count(kh.kmer_hash, seed_val));
		Vector<Long> bifurcation = new Vector<Long>();
		// String right = forward_extend(kh.list.get(0).getKey(), bifurcation,
		// kh);
		String right = forward_extend(seed_val, bifurcation, kh);
		String left = reverse_extend(seed_val, bifurcation, kh);
		String trunk = left + right.substring(kh.kmer_length);
		Node trunkN = new Node();
		trunkN.sequence = trunk;
		if (trunkN.sequence.length() < 200) {
			return false;
		}
		SplicingGraph splicing_graph = new SplicingGraph();
		int p = splicing_graph.add_node(trunkN);
		// System.out.println("******************" + p);
		// // System.out.println("原序列："+splicing_graph.node_set.size());
		// System.out.println("原序列" + " size:" +
		// splicing_graph.node_set.get(0).sequence.length() + " "
		// + splicing_graph.node_set.get(0).sequence);
		splicing_graph.extend_again(kh, bifurcation);
		// System.out.println("第一次扩展:" + " size:" +
		// splicing_graph.node_set.get(0).sequence.length() + " "
		// + splicing_graph.node_set.get(0).sequence);
		// System.out.println();
		splicing_graph.forward_extend_use_pairInfo(kh, load_Read.read_vector, p, bifurcation);
		// System.out.println("paired——end扩展：" + " size:" +
		// splicing_graph.node_set.get(0).sequence.length() + " "
		// + splicing_graph.node_set.get(0).sequence);
		splicing_graph.reverse_extend_use_pairInfo(kh, load_Read.read_vector, p, bifurcation);

		node_jihe.add(splicing_graph.node_set.get(p).sequence);

		System.out.println();
		System.out.println(bifurcation);
		return true;
	}

	public void rewrite_nodeSet(Set node_jihe) {
		// TODO Auto-generated method stub
		node_set.clear();
		Iterator<String> it = node_jihe.iterator();
		while (it.hasNext()) {
			String node = it.next();
			Node node1 = new Node();
			node1.sequence = node;
			node_set.add(node1);
		}
	}

	// public boolean bulid_graph(kmerHash kh, long seed_val) {
	//
	// if (!init_trunk(kh, seed_val)) {
	// return false;
	// }
	//
	// return false;
	//
	// }

	public void forward_check_and_extend(kmerHash kh, int node_index) {
		int length = node_set.get(node_index).sequence.length() - kh.kmer_length;
		Vector<Long> bifurcations = new Vector<Long>();
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		if (length < 0) {
			return;
		}
		for (int i = 0; i < length; i++) {
			String kmer = node_set.get(node_index).sequence.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
			if (candidates.size() == 0) {
				continue;
			}
			long candidate = 0l;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					bifurcations.add(intval);
					break;
				}
			}

		}
		grow_and_branch(kmer_map, data, p, bifurcation_points);
	}

	public void grow_and_branch(kmerHash kh, int node_index, Vector<Long> bifurcations) {
		while (bifurcations.size() > 0) {
			long intval = bifurcations.lastElement();
			bifurcations.remove(bifurcations.size() - 1);
			Vector<Long> bifurcations_more = new Vector<Long>();
			String str = forward_extend(intval, bifurcations_more, kh);
			if (bifurcations_more.size() > 0 && bifurcations_more.get(0) == intval) {
				// 如果intval还有两个或以上没用的candidate
				bifurcations.add(intval);
				bifurcations_more.remove(0);
			}
			String end_kmer = str.substring(str.length() - kh.kmer_length); // str的最后一个kmer
			long end_intval = baseOptions.kmerToIntval(end_kmer);
			List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
			candidates = kh.get_forward_candidates(end_intval, kh.kmer_hash);
			if (candidates.size() > 0) {
				int bubble_val = add_bubble(node_index, str, kh);
				if(bubble_val>0){
					forward_branches.add(bubble_val);
				}
				if(bubble_val==-2){
					
				}
			}

			if (candidates.size() > 0) { // add bubble

				node_idx_t r = add_bubble(p, sequence, kmer_map);
				if (r > 0) {
					forward_branches.insert(r);
				} else if (r == -2) { // return -2 if branch //找不到q
					r = add_branch(p, sequence, kmer_map, data);
					if (r > 0) {
						forward_branches.insert(r);
						node_set_[r].sequence = node_set_[r].sequence.substr(g_kmer_length);
					}
				}
			} else { // add a branch

				node_idx_t r = add_branch(p, sequence, kmer_map, data);
				if (r > 0) {
					forward_branches.insert(r);
					// change sequence of node r
					node_set_[r].sequence = node_set_[r].sequence.substr(g_kmer_length);
				}
			}

		} // while
	}

	public int add_bubble(int node_p, String str, kmerHash kh) {

	    if (str.length() < 2 * kh.min_anchor_length)
	    {
	    	return -1;
	    }
		String anchor_left=str.substring(0,kh.kmer_length);
		int start=node_set.get(node_p).sequence.indexOf(anchor_left);
		//如果找不到则不断缩小anchor_left
		if(start==-1){
			while(anchor_left.length()>kh.min_anchor_length){
				anchor_left=anchor_left.substring(0,anchor_left.length()-1);
				start=node_set.get(node_p).sequence.indexOf(anchor_left);
				if(start!=-1){
					break;    //直至找到为止
				}
			}
			if(start==-1){
				System.out.println("无法找到起始anchor！");
				return -1;
			}
		}
		
	    String anchor_right=str.substring(str.length()-kh.kmer_length+1);
	    int node_q=-1;
	    while(anchor_right.length()>kh.min_anchor_length){
	    	node_q=find_node_index(anchor_right);
	    	if(node_q>=0){
	    		break;
	    	}
	    	anchor_right=anchor_right.substring(1);
	    }
	    if(node_q==-1){
	    	restore_kmers(str, kh);
	    	return -2;
	    }
		
	int anchor_left_length = anchor_left.length(); 
	int anchor_right_length = anchor_right.length();
    int end=node_set.get(node_q).sequence.indexOf(anchor_right);

    //根据连接处的cov判断是否可以连接
    String jun1=str.substring((int)(0.5*anchor_left_length),(int)(0.5*anchor_left_length+kh.kmer_length));
    long jun1_cov=kh.get_readset_count(kh.kmer_hash, baseOptions.kmerToIntval(jun1));
    String jun2=str.substring(str.length()-kh.kmer_length-(int)(0.5*anchor_right_length));
    long jun2_cov=kh.get_readset_count(kh.kmer_hash, baseOptions.kmerToIntval(jun2));
    if(jun1_cov<kh.min_function_count||jun2_cov<kh.min_function_count){
    	System.out.println("连接点cov太小，无法连接！");
    	return -1;
    }
	
    int length=str.length()-anchor_left_length-anchor_right_length;   //str有多少未被使用

    if(node_p==node_q){   //顶点本身回归
    	// TODO Auto-generated method stub    可修改
    	if(end<=start){
    		return -1;
    	}
    	int distance=end-start-anchor_left_length;  //表示end与start之间是否比anchor_left大
	    if(length+distance<4){
	    	return -1;
	    }
	    if ((distance == 0 && length < kh.kmer_length) || (length == 0 && distance < kh.kmer_length))
	    {
	    	return -1;
	    }
	    if ( distance > 0 && length == distance
		        && is_similar(node_set.get(node_p).sequence.substring(start+anchor_left_length, start+anchor_left_length+distance),
		        str.substring(anchor_left_length,anchor_left_length+length), 'F')
		        ) {
					//重合
		          return -1;
		      }
	    if ( length <= 0 && distance <= 0) {
			return -1;
		      } else if ( length < 0) {
			anchor_left_length = anchor_left_length + length;
		      } else if (distance < 0) {
			anchor_left_length = anchor_left_length + distance;  //设置 anchor_left与anchor_right刚好完全占据str
		      }
	    if ((int)start+anchor_left_length <= kh.kmer_length) // necessary 
			{
	    	return -1;
			}
	    Node node1=new Node();
	    Node node2=new Node();
	    node1.sequence=node_set.get(node_p).sequence.substring(end);
        int node1_index=add_node(node1);
        node_set.get(node1_index).children=node_set.get(node_p).children;
        node_set.get(node_p).children.clear();
        int node2_index=-1;
        if(distance<=0){ //此时length>0
        	node_set.get(node_p).sequence=node_set.get(node_p).sequence.substring(0,start+anchor_left_length);
        	node_set.get(node_p).addChildren(node1_index);
        }
        else{
        	Node node3=new Node();
        	if(length<0){
        		node3.sequence=node_set.get(node_p).sequence.substring(start+anchor_left_length,start+anchor_left_length+distance-length);
        	}else{
        		node3.sequence=node_set.get(node_p).sequence.substring(start+anchor_left_length,start+anchor_left_length+distance);
        	}
        	int node3_index=add_node(node3);
        	node_set.get(node3_index).addChildren(node1_index);
        	node_set.get(node_p).sequence.substring(0,start+anchor_left_length);
        	node_set.get(node_p).addChildren(node3_index);
        }
        if(length>0){
        	if(distance<0){
        		node2.sequence=str.substring(anchor_left_length,anchor_left_length+length-distance);
        	}
        	else{
        		node2.sequence=str.substring(anchor_left_length,anchor_left_length+length);
        	}
        	node2_index=add_node(node2);
        	node_set.get(node2_index).addChildren(node1_index);
        	node_set.get(node_p).addChildren(node2_index);
        }
        else{
        	node_set.get(node_p).addChildren(node1_index);
        }
        return node2_index;   //新增顶点标号
    }
    else{
    	Set<Integer> checked=new HashSet<Integer>();
    	if(node_q==0||has_path(node_p,node_q,checked)){
    		return -1;
    	}
    	if(node_set.get(node_p).isChild(node_q)&&length>0&&node_set.get(node_p).sequence.length()-start+end<kh.kmer_length){
    		return -1;
    	}
    	Node node1=new Node();
    	Node node2=new Node();
    	int node1_index=-1;
    	int node2_index=-1;
    	if(length<0){
    		anchor_left_length = anchor_left_length + length;
    	}
    	if(start+anchor_left_length<=kh.kmer_length){
    		return -1;
    	}
    	if(start+anchor_left_length<node_set.get(node_p).sequence.length()){
    		//只能小于或等于  如果等于，无需再进行操作
    		node1.sequence=node_set.get(node_p).sequence.substring(start+anchor_left_length);
    		node_set.get(node_p).sequence.substring(0,start+anchor_left_length);
    		node1_index=add_node(node1);
    		node_set.get(node1_index).children=node_set.get(node_p).children;
    		node_set.get(node_p).children.clear();
    		node_set.get(node_p).addChildren(node1_index);
    	}
    	if(length>0){
    		node2.sequence=str.substring(anchor_left_length,length);
    		node2_index=add_node(node2);
    		node_set.get(node_p).addChildren(node2_index);
    	}else{
    		node2_index=node_p;  //?
    	}
    	if(end==0){
    		node_set.get(node2_index).addChildren(node_q);
    	}
    	else{
    		Node node3 =new Node();
    		if(node_set.get(node_q).sequence.length()>=end){
    			node3.sequence=node_set.get(node_q).sequence.substring(end);
    			int node3_index=add_node(node3);
    			node_set.get(node3_index).children=node_set.get(node_q).children;
    			node_set.get(node_q).children.clear();
    			node_set.get(node_q).addChildren(node3_index);
    			node_set.get(node2_index).addChildren(node3_index);
    			node_set.get(node_q).sequence=node_set.get(node_q).sequence.substring(0,end);
    		}
    	}
    	if(node2_index==node_p){
    		return -1;
    	}
    	else{
    		return node2_index;
    	}
    }
	
}

	public boolean has_path(int node_p, int node_q, Set checked) {
		// 判断q到p是否有路径
		boolean flag = false;
		if (node_set.get(node_q).children.size() == 0 || checked.contains(node_q)) {
			return false;
		} else {
			checked.add(node_q);
		}
		for (int i = 0; i < node_set.get(node_q).children.size(); i++) {
			if (node_set.get(node_q).children.get(i) == node_p) {
				flag = true;
				break;
			} else {
				flag = has_path(node_p, node_set.get(node_q).children.get(i), checked);
				if (flag) {
					break;
				}
			}
		}
		return flag;
	}

	public int find_node_index(String anchor) {
		int idx = -1;
		for (int i = 0; i < node_set.size(); ++i) {
			if (node_set.get(i).sequence.indexOf(anchor) != -1) {
				idx = i;
				break;
			}
		}
		return idx;
	}

	public void restore_kmers(String str, kmerHash kh) {

		if (str.length() < kh.kmer_length) {
			return;
		}

		for (int i = 0; i <= str.length() - kh.kmer_length; i++) {
			String kmer = str.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			restore_kmers.add(intval);
		}

	}

	
	//如果在扩展过程中找不到对应的q  则回不到主干，只能作为分支
	  public int add_branch(int node_p,String branch,kmerHash kh) {

	    if (node_set.get(node_p).sequence.length() <= kh.kmer_length)
	    {
	    	return -1;
	    }
		if(branch.length()>kh.kmer_length+kh.min_exon_length){
			String start_kmer=branch.substring(0,kh.kmer_length);
			int start=node_set.get(node_p).sequence.indexOf(start_kmer);
			if(start==-1){
				restore_kmers(branch.substring(1), kh);
				return -1;
			}
			if(node_p==0&&start<kh.pair_gap_length){
				restore_kmers(branch.substring(1), kh);
				return -1;
			}
			if (kh.is_paired_end) { 

		        long count = kh.get_readset_count(kh.kmer_hash, baseOptions.kmerToIntval(start_kmer));
		        boolean is_branch = check_forward_branch_with_pair_info(kh, branch, count);
		        
		        if (!is_branch) {
		          restore_kmers(branch.substring(1), kh);
		          return -1;
		        }
		      }
			Node node1=new Node();
			Node node2=new Node();
			if(start+kh.kmer_length<node_set.get(node_p).sequence.length()){
				//需要拆分
				node1.sequence=node_set.get(node_p).sequence.substring(start+kh.kmer_length);
			}
			node2.sequence=branch;
			if(is_similar(branch.substring(kh.kmer_length),node1.sequence,'F')){
				if(node_set.get(node_p).children.size()==0&&node1.sequence.length()<kh.min_exon_length){
					return -2;
				}
				else if(node1.sequence.length()+kh.kmer_length>branch.length()){
					return -1;
				}
			}
			if(node_set.size()==1&&node1.sequence.length()<2*kh.kmer_length){
				node_set.get(node_p).sequence=node_set.get(node_p).sequence.substring(0,start)+branch;
				return -2;
			}
			if(node1.sequence.length()!=0){
				int node1_index=add_node(node1);
				node_set.get(node1_index).children=node_set.get(node_p).children;
				node_set.get(node_p).children.clear();
				node_set.get(node_p).addChildren(node1_index);
			}
			int node2_index=add_node(node2);
			node_set.get(node_p).addChildren(node2_index);
			node_set.get(node_p).sequence=node_set.get(node_p).sequence.substring(0, start+kh.kmer_length);
			return node2_index;
		}
		else{
			return -1;
		}
	    
	  }

	  public boolean check_forward_branch_with_pair_info(kmerHash kh,String branch, long count) {

		    boolean is_branch = false;
		    String check=branch.substring(kh.kmer_length,3*kh.kmer_length);
		   int middle_read_id=load_Read.read_vector.size()/2;
		    Set<Integer> reads= new HashSet<Integer>();
		    set_reads(kh, check, reads);
		    Set<Long> kmers_in_branch =new HashSet<Long>();
		    int support = 0;
		    int unsupport = 0;
		    Iterator<Integer> it = reads.iterator();
			int paired_end_read_id = 0;
			String extend_str = "";
			while (it.hasNext()) {
				int read_id = it.next();
				if (kh.fr_strand == 1) { // 2-> <-1 此时应该把1反过来 根据2找1
					if (read_id < middle_read_id ) {
						String mate_left=load_Read.read_vector.get(read_id);
						if (compatible(branch, mate_left, kh)) {
				              String mate_right = load_Read.read_vector.get(read_id+middle_read_id);
				              int used = 0;
				              for (int i = 0; i <= mate_right.length()-kh.kmer_length; i++) {
				                long intval=baseOptions.kmerToIntval(mate_right.substring(i,i+kh.kmer_length));
				                if (has_been_used(intval))
				                  used++;
				              }
				              if (used == 0) {
				                unsupport++;
				                if ((unsupport >= kh.min_ratio_welds * count) && (unsupport >= kh.min_reads_span_junction))
				                  break;
				              } else if (used + 5 >= mate_right.length() - kh.kmer_length) {
				                support++;
				                if ((support >= kh.min_ratio_welds * count) && (support >= kh.min_reads_span_junction)) {
				                  is_branch = true;
				                  break;
				                }
				              }
				            }
					} 
				} else if (kh.fr_strand == 2) {//->1    <-2
					if (read_id >= middle_read_id) {  
			           String mate_right = load_Read.read_vector.get(read_id);
			            if (compatible(branch, mate_right, kh)) {
			              String mate_left =load_Read.read_vector.get(read_id-middle_read_id);
			              int used = 0;
			              for (int i = 0; i <= (int)mate_left.length() - kh.kmer_length; ++i) {
			                long intval = baseOptions.kmerToIntval(mate_left.substring(i,kh.kmer_length));
			                if (has_been_used(intval))
			                  used++;
			              }
			              if (used == 0) {
			                unsupport++;
			                if ((unsupport >= kh.min_ratio_welds * count) && (unsupport >= kh.min_reads_span_junction))
			                  break;
			              } else if (used + 5 >= (int)mate_left.length() - kh.kmer_length) {
			                support++;
			                if ((support >= kh.min_ratio_welds * count) && (support >= kh.min_reads_span_junction)) {
			                  is_branch = true;
			                  break;
			                }
			              }
			            }
			          }				
					}

			   } // else
	   

	    return is_branch;
	  }
		    
}
