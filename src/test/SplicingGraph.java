package test;

import java.util.List;
import java.util.Vector;

public class SplicingGraph {
	public class Node {
		private Vector<Long> parents = new Vector<Long>();
		private Vector<Long> children = new Vector<Long>();
		private String sequence;

		public boolean addParents(long parent) {
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

		public boolean addChildren(long child) {
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

		public boolean isChild(long child) {
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

		public boolean isParent(long parent) {
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

	public void add_node() {
		Node node = new Node();
		node.setSequence("zhangjianxi");
		System.out.println(node.getSequence());
	}

}
