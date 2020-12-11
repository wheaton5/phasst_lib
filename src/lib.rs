extern crate hashbrown;
extern crate flate2;
extern crate byteorder;
use flate2::read::GzDecoder;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
use std::fs::File;
use std::io::{BufWriter, Write};

use byteorder::{ByteOrder, LittleEndian};

use hashbrown::{HashMap, HashSet};
use std::str;





fn eat_i32<R: Read>(reader: &mut BufReader<R>, buf: &mut [u8;4]) -> Option<i32> {
    let bytes = reader.read(buf).expect("could not read binary molecules file");
    if bytes < 4 { return None; }
    Some(LittleEndian::read_i32(buf))
}

#[allow(dead_code)]
pub struct Molecules {
    //long_read_molecule_list: Vec<i32>,
    longread_molid_offsets: Vec<i32>,
    linked_read_molecules: HashMap<i32, Vec<i32>>,
    hic_molecules: HashMap<i32, Vec<i32>>,
    long_read_molecules: HashMap<i32, Vec<i32>>,
    long_read_molecule_list: HashMap<i32, Vec<i32>>,
    long_read_het_molecules: HashMap<i32, Vec<i32>>,
    long_read_het_molecule_list: HashMap<i32, Vec<i32>>,
    long_read_molecule_positions: HashMap<i32, Vec<i32>>,
    hom_linked_read_molecules: HashMap<i32, HashSet<i32>>,
    hom_long_read_molecules: HashMap<i32, HashSet<i32>>,
    hom_hic_molecules: HashMap<i32, HashSet<i32>>,
    het_linked_read_molecules: HashMap<i32, HashSet<i32>>,
    het_long_read_molecules: HashMap<i32, HashSet<i32>>,
    het_hic_molecules: HashMap<i32, HashSet<i32>>,
}

impl Molecules {
    pub fn get_linked_read_variants<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.linked_read_molecules.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_long_read_molecules<'a>(&'a self) -> Box<dyn Iterator<Item=&i32>+'a> {
        Box::new(self.long_read_molecules.keys())
    }
    pub fn get_linked_read_molecules<'a>(&'a self) -> Box<dyn Iterator<Item=&i32>+'a> {
        Box::new(self.linked_read_molecules.keys())
    }

    pub fn get_hic_molecules<'a>(&'a self) -> Box<dyn Iterator<Item=&i32>+'a> {
        Box::new(self.hic_molecules.keys())
    }
    pub fn get_long_read_variants<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.long_read_molecules.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_long_read_variants_ordered<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.long_read_molecule_list.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_long_read_variants_and_positions<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=(&i32, &i32)>+'a> {
        match self.long_read_molecule_list.get(&mol.abs()) {

            Some(x) => {
                match self.long_read_molecule_positions.get(&mol.abs()) {
                    Some(y) => Box::new(x.iter().zip(y.iter())),
                    None => Box::new(std::iter::empty()),
                }
            },
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_hic_variants<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hic_molecules.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_hom_linked_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_linked_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_hom_long_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_long_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    pub fn get_hom_hic_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_hic_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    pub fn get_het_linked_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_linked_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_het_long_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_long_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_het_hic_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_long_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_variants<'a>(&'a self, i: &i32, kmer_type: KmerType) -> Box<dyn Iterator<Item=&i32>+'a> {
        match kmer_type {
            KmerType::PairedHet => Box::new(
                    self.get_linked_read_variants(*i).chain(
                    self.get_long_read_variants(*i)).chain(
                    self.get_hic_variants(*i))),
            KmerType::UnpairedHet => Box::new(
                    self.get_het_linked_read_variants(i).chain(
                    self.get_het_long_read_variants(i)).chain(
                    self.get_het_hic_variants(i))),
            KmerType::Homozygous => Box::new(
                    self.get_hom_linked_read_variants(i).chain(
                    self.get_hom_long_read_variants(i)).chain(
                    self.get_hom_hic_variants(i))),
        }
    }
}

impl Variants {
    pub fn get_variant_iter<'a>(&'a self, kmer_type: KmerType) -> Box<dyn Iterator<Item=&i32>+'a> {
        match kmer_type {
            KmerType::PairedHet => Box::new(self.paired_het_variants.iter()),
            KmerType::UnpairedHet => Box::new(std::iter::empty()),
            KmerType::Homozygous => Box::new(std::iter::empty()),
        }
    }
    pub fn len(&self) -> usize {
        return self.linked_read_variants.len();
    }
    pub fn get_linked_read_molecules(&self, i: i32) -> std::slice::Iter<'_, i32> {
        return self.linked_read_variants[(i.abs() - 1) as usize].iter();
    }
    pub fn get_num_linked_read_molecules(&self, i: i32) -> usize {
        return self.linked_read_variants[(i.abs() - 1) as usize].len();
    }
    pub fn get_long_read_molecules(&self, i: i32) -> std::slice::Iter<'_, i32> {
        return self.long_read_variants[(i.abs() - 1) as usize].iter();
    }
    pub fn get_num_long_read_molecules(&self, i: i32) -> usize {
        return self.long_read_variants[(i.abs() - 1) as usize].len();
    }
    pub fn get_hic_molecules(&self, i: i32) -> std::slice::Iter<'_, i32> {
        return self.hic_variants[(i.abs() - 1) as usize].iter();
    }
    pub fn get_num_hic_molecules(&self, i: i32) -> usize {
        return self.hic_variants[(i.abs() - 1) as usize].len();
    }
    #[allow(dead_code)]
    pub fn get_hom_linked_read_molecules<'a>(&'a self, i: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_linked_read_variants.get(&i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    #[allow(dead_code)]
    pub fn get_hom_long_read_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_long_read_variants.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    #[allow(dead_code)]
    pub fn get_het_linked_read_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_linked_read_variants.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    #[allow(dead_code)]
    pub fn get_het_long_read_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_long_read_variants.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    pub fn get_molecules<'a>(&'a self, i: &i32, data_type: DataType) 
            -> Box<dyn Iterator<Item=&i32>+'a> {
        match data_type {
            DataType::HC => Box::new(self.get_hic_molecules(*i)),
            DataType::PB => Box::new(self.get_long_read_molecules(*i)),
            DataType::LR => Box::new(self.get_linked_read_molecules(*i)),
            DataType::PBLR => Box::new(self.get_long_read_molecules(*i).chain(self.get_linked_read_molecules(*i))),
        }
    }
    pub fn get_num_molecules(&self, i: &i32, data_type: DataType) -> usize {
        match data_type {
            DataType::HC => self.get_num_hic_molecules(*i),
            DataType::PB => self.get_num_long_read_molecules(*i),
            DataType::LR => self.get_num_linked_read_molecules(*i),
            DataType::PBLR => self.get_num_long_read_molecules(*i) + self.get_num_linked_read_molecules(*i),
        }
    }
}

#[allow(dead_code)]
pub enum DataType {
    HC, // hic
    PB, // pacbio
    LR, // linked reads
    PBLR, //pb and lr
}

struct Crib {
    variants: HashMap<i32, (i32, usize, usize)>, // map from kmer_id to crib molecule id and number seen and position
    molecules: HashMap<i32, HashMap<i32, (usize, usize)>>, // map from crib molecule id to a map from kmer_id to order index
}

fn load_crib(crib: &String, kmers: &Kmers) -> Crib {
    let mut mol_id = 1;
    let mut variants: HashMap<i32, (i32, usize, usize)> = HashMap::new();
    let mut molecules: HashMap<i32, HashMap<i32, (usize, usize)>> = HashMap::new();

    let mut paired_het_variants: HashSet<i32> = HashSet::new();
    let mut long_read_variants: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecule_list: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut hom_long_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut bufi32 = [0u8; 4];
    let mut buf2 = [0u8; 4];
    let f = File::open(crib)
        .expect(&format!("Unable to open crib {}", crib));
    let mut reader = BufReader::new(f);
    let mut ordered_kmers: Vec<(i32, i32, usize)> = Vec::new();
    'outer: loop { // ok here we go again. Pacbio/longread data. Format is i32s until you hit a zero. when you hit two zeros you are done
        let mut vars: HashSet<i32> = HashSet::new();
        let mut varlist: Vec<i32> = Vec::new();
        let mut hom_kmers: HashSet<i32> = HashSet::new();
        //let mut het_kmers: HashSet<i32> = HashSet::new();
        let mut order = 0;
        loop {
            if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                if kmer_id == 0 { break; }
                if let Some(position) = eat_i32(&mut reader, &mut buf2){
                    let position = position + 10;
                    //eprintln!("position {}", position);
                    let mut molecule = 0;
                    let mut number = 0;
                    let mut position1 = 0;
                    let mut has = false;
                   

                    
                    let var_order = molecules.entry(mol_id).or_insert(HashMap::new());
                    var_order.insert(kmer_id.abs(), (order, position as usize));
                    match kmers.kmer_type.get(&kmer_id) {
                        Some(KmerType::PairedHet) => {
                            ordered_kmers.push((mol_id, kmer_id.abs(), position as usize));
                            //eprintln!("paired het");
                            vars.insert(kmer_id);
                            paired_het_variants.insert(kmer_id.abs());
                            varlist.push(kmer_id);
                             if let Some((mol, num, pos)) = variants.get(&kmer_id.abs()) {
                                molecule = *mol;
                                number = *num;
                                position1 = *pos;
                                has = true;
                            } else {
                                variants.insert(kmer_id.abs(), (mol_id, 1, position as usize));
                            }
                            if has {
                                variants.insert(kmer_id.abs(), (mol_id, number+1, position as usize));
                            }
                        },
                        Some(KmerType::UnpairedHet) => {
                            //het_kmers.insert(kmer_id);
                        },
                        Some(KmerType::Homozygous) => {
                            hom_kmers.insert(kmer_id);
                        },
                        None => { eprintln!("no kmer type? {}", kmer_id); }
                    }
                }
            } else { break 'outer; }
            order += 1;
        }
            for kmer_id in vars.iter() {
                let var_bcs = long_read_variants.entry(kmer_id.abs()).or_insert(HashSet::new());
                if kmer_id < &0 { var_bcs.insert(-mol_id); } else { var_bcs.insert(mol_id); }
            }
            long_read_molecules.insert(mol_id, vars); 
            long_read_molecule_list.insert(mol_id, varlist);
            
        if hom_kmers.len() > 0 {
            for kmer_id in hom_kmers.iter() {
                let var_bcs = hom_long_read_kmers.entry(*kmer_id).or_insert(HashSet::new());
                var_bcs.insert(mol_id);
            }
        }

        hom_long_read_molecules.insert(mol_id, hom_kmers);
        mol_id += 1;
        
    }

    Crib{
        variants: variants,
        molecules: molecules,
    }
}

#[allow(dead_code)]
pub struct Variants {
    paired_het_variants: HashSet<i32>,
    linked_read_variants: Vec<Vec<i32>>,
    hic_variants: Vec<Vec<i32>>,
    long_read_variants: Vec<Vec<i32>>,
    long_read_het_variants: Vec<Vec<i32>>,
    hom_linked_read_variants: HashMap<i32, HashSet<i32>>,
    hom_long_read_variants: HashMap<i32, HashSet<i32>>,
    hom_hic_variants: HashMap<i32, HashSet<i32>>,
    het_linked_read_variants: HashMap<i32, HashSet<i32>>,
    het_long_read_variants: HashMap<i32, HashSet<i32>>,
    het_hic_variants: HashMap<i32, HashSet<i32>>,
}

pub fn load_molecule_kmers(txg_mols: &Vec<String>, hic_mols: &Vec<String>, longread_mols: &Vec<String>, kmers: &Kmers) -> (Variants, Molecules){
    let mut linked_read_variants: HashMap<i32, HashSet<i32>> = HashMap::new(); // map from variant_id to list of molecule_ids
    let _hic_variants: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_variants: HashMap<i32, HashSet<i32>> = HashMap::new();
    
    let mut linked_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new(); //map from molecule_id to list of variant_ids
    let mut hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecule_list: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut long_read_molecule_positions: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut long_read_het_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_het_molecule_list: HashMap<i32, Vec<i32>> = HashMap::new();

    let mut paired_het_variants: HashSet<i32> = HashSet::new();

    let mut hom_linked_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_long_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_linked_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let hom_hic_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let hom_hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();

    let het_linked_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_long_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_linked_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_hic_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();

    let mut bufi32 = [0u8; 4];
    let mut buf2 = [0u8; 4];
    let mut max_var = 0;
    let mut max_molid = 0;
    for txg_file in txg_mols.iter() {
        println!("{}",txg_file);
        let f = File::open(txg_file.to_string())
            .expect(&format!("Unable to open txg file {}", txg_file));
        let mut reader = BufReader::new(f);
        loop {
            if let Some(barcode_id) = eat_i32(&mut reader, &mut bufi32) {
                if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                    if kmer_id.abs() > max_var { max_var = kmer_id.abs() }
                    if barcode_id > max_molid { max_molid = barcode_id; }
                    match kmers.kmer_type.get(&kmer_id).unwrap() {
                        KmerType::PairedHet => {
                            let bc_vars = linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            bc_vars.insert(kmer_id);
                            paired_het_variants.insert(kmer_id.abs());
                        },
                        KmerType::UnpairedHet => {
                            //eprintln!("unpaired het kmer");
                            //let bc_vars = het_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            //bc_vars.insert(kmer_id);
                            ()
                        },
                        KmerType::Homozygous => {
                            let bc_vars = hom_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            bc_vars.insert(kmer_id);
                            ()
                        },
                    }
                    
                } else { break; }
            } else { break; }
        }
    }
    
    eprintln!("reducing to good linked read molecules with > 5 variants, right now we have {}", linked_read_molecules.len());
    //linked_read_molecules.retain(|_key, value| value.len() > 10 && value.len() < 5000); 
    for (mol, varset) in linked_read_molecules.iter() {
        for var in varset.iter() {
            let var_bcs = linked_read_variants.entry(var.abs()).or_insert(HashSet::new());
            if *var < 0 { var_bcs.insert(-mol); } else { var_bcs.insert(*mol); }
        }
    }
    
    eprintln!("{} good linked read molecules", linked_read_molecules.len());
   
    for (mol, varset) in hom_linked_read_molecules.iter() {
        for var in varset.iter() {
            let var_bcs = hom_linked_read_kmers.entry(*var).or_insert(HashSet::new());
            var_bcs.insert(*mol);
        }
    }
    
    eprintln!("{} hom linked read molecules, {} hom linked read kmers",hom_linked_read_molecules.len(), hom_linked_read_kmers.len());
    //molecules.clear();

    /*
    let mut bad_vars = 0;
    let mut vars_to_remove: Vec<i32> = Vec::new();
    let mut bad_var_set: HashSet<i32> = HashSet::new();
    for (var, var_mols) in linked_read_variants.iter() {
        let mut mols: HashSet<i32> = HashSet::new();
        for mol in var_mols.iter() {
            mols.insert(mol.abs());
        }
        if mols.len() + 1 < var_mols.len() {
            bad_vars += 1;
            //println!("bad variant {}, {} mols with both versions, {} total mols with either.", kmer_id_to_kmer.get(var).unwrap(), var_mols.len() - mols.len(), mols.len());
            for mol in var_mols.iter() {
                if let Some(linked_read_molecule) = linked_read_molecules.get_mut(&mol.abs()) {
                    linked_read_molecules.remove(var);
                    linked_read_molecules.remove(&-var);
                }
            }
            vars_to_remove.push(var.abs());
            bad_var_set.insert(var.abs());
        }
    }
    */
    let mut mol_id = max_molid + 1;
    for hic_file in hic_mols.iter() {
        let f = File::open(hic_file.to_string())
            .expect(&format!("Unable to open hic file {}", hic_file));
        let mut reader = BufReader::new(f);
        'outerhic: loop { // now deal with hic data, format is i32s until you hit a 0 if i get to 2 0's we are done
            //break 'outerhic;
            let mut vars: HashSet<i32> = HashSet::new();
            let mut vars2: Vec<i32> = Vec::new();
            loop {
                if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                    //println!("kmer id {}", kmer_id);

                    if kmer_id == 0 { if vars.len() == 0 { break 'outerhic; } else { break; } }
                    match kmers.kmer_type.get(&kmer_id).unwrap() {
                        KmerType::PairedHet => {
                            //eprintln!("paired het kmer");
                                vars.insert(kmer_id); 
                                vars2.push(kmer_id);
                            ()
                        },
                        KmerType::UnpairedHet => {
                            //eprintln!("unpaired het kmer");
                            //vars.insert(kmer_id); 
                            vars2.push(kmer_id);
                            //let bc_vars = het_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            //bc_vars.insert(kmer_id);
                            ()
                        },
                        KmerType::Homozygous => {
                            //eprintln!("homozygous kmer");
                            //vars.insert(kmer_id); 
                            vars2.push(kmer_id);
                            //let bc_vars = hom_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            //bc_vars.insert(kmer_id);
                            ()
                        },
                    }
                    //if !bad_var_set.contains(&kmer_id.abs()) {
                    
                    
                    if kmer_id.abs() > max_var { max_var = kmer_id.abs() }
                    //}
                } else { break 'outerhic; }
            }
            mol_id += 1;
            if vars.len() > 1 {
                hic_molecules.insert(mol_id, vars);
            }
        }
    }

    eprintln!("num hic molecules is {}", hic_molecules.len());
    
    mol_id += 1; 
    let mut longread_mol_id_starts: Vec<i32> = Vec::new();
    for longread_file in longread_mols.iter() {
        longread_mol_id_starts.push(mol_id);
        let f = File::open(longread_file.to_string())
            .expect(&format!("Unable to open longread file {}", longread_file));
        let mut reader = BufReader::new(f);
        'outer: loop { // ok here we go again. Pacbio/longread data. Format is i32s until you hit a zero. when you hit two zeros you are done
            let mut vars: HashSet<i32> = HashSet::new();
            let mut varlist: Vec<i32> = Vec::new();
            let mut varlist_positions: Vec<i32> = Vec::new();
            let mut hom_kmers: HashSet<i32> = HashSet::new();
            let mut het_kmers: HashSet<i32> = HashSet::new();
            let mut het_kmers_list: Vec<i32> = Vec::new();
            loop {
                if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                    
                    if kmer_id == 0 { break; }
                    if let Some(position) = eat_i32(&mut reader, &mut buf2) {
                        let position = position + 10;
                        if kmer_id.abs() > max_var { max_var = kmer_id.abs(); }
                        match kmers.kmer_type.get(&kmer_id) {
                            Some(KmerType::PairedHet) => {
                                vars.insert(kmer_id);
                                paired_het_variants.insert(kmer_id.abs());
                                varlist.push(kmer_id);
                                varlist_positions.push(position);
                            },
                            Some(KmerType::UnpairedHet) => {
                                het_kmers.insert(kmer_id);
                                het_kmers_list.push(kmer_id);
                            },
                            Some(KmerType::Homozygous) => {
                                hom_kmers.insert(kmer_id);
                            },
                            None => { eprintln!("no kmer type? {}", kmer_id); }
                        }
                    }
                } else { break 'outer; }
            }
                for kmer_id in vars.iter() {
                    let var_bcs = long_read_variants.entry(kmer_id.abs()).or_insert(HashSet::new());
                    if kmer_id < &0 { var_bcs.insert(-mol_id); } else { var_bcs.insert(mol_id); }
                }
                long_read_molecules.insert(mol_id, vars); 
                long_read_molecule_list.insert(mol_id, varlist);
                long_read_molecule_positions.insert(mol_id, varlist_positions);
                long_read_het_molecules.insert(mol_id, het_kmers);
                long_read_het_molecule_list.insert(mol_id, het_kmers_list);
                
                
            /*
            if het_kmers.len() > 0 {
                for kmer_id in het_kmers.iter() {
                    let var_bcs = het_long_read_kmers.entry(*kmer_id).or_insert(HashSet::new());
                    var_bcs.insert(mol_id);
                }
            }
            */
            if hom_kmers.len() > 0 {
                for kmer_id in hom_kmers.iter() {
                    let var_bcs = hom_long_read_kmers.entry(*kmer_id).or_insert(HashSet::new());
                    var_bcs.insert(mol_id);
                }
            }
            hom_long_read_molecules.insert(mol_id, hom_kmers);
            mol_id += 1;
            
        }
    }
    
    //long_read_molecules.retain(|_key, value| value.len() > 1 && value.len() < 5000); 
    eprintln!("{} good long read molecules", long_read_molecules.len());
    
    eprintln!("{} hom long read kmers, {} hom read long molecules", hom_long_read_kmers.len(), hom_long_read_molecules.len());

    

    let mut txg_vars: Vec<Vec<i32>> = Vec::new();
    let mut hic_vars: Vec<Vec<i32>> = Vec::new();
    let mut pb_vars: Vec<Vec<i32>> = Vec::new();
    let mut pb_het_vars: Vec<Vec<i32>> = Vec::new();
    for _i in 0..(max_var+1) {
        txg_vars.push(Vec::new());
        hic_vars.push(Vec::new());
        pb_vars.push(Vec::new());
        pb_het_vars.push(Vec::new());
    }

    let mut txg_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut hic_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut pb_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut pb_het_mols: HashMap<i32, Vec<i32>> = HashMap::new();

    for (mol_id, kmer_ids) in linked_read_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                txg_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { txg_vars[(var - 1) as usize].push(mol_id) }
        }

        txg_mols.insert(mol_id, kmer_ids_sorted);
    }

    for (mol_id, kmer_ids) in hic_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                hic_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { hic_vars[(var - 1) as usize].push(mol_id) }
        }
        hic_mols.insert(mol_id, kmer_ids_sorted);
    }
    for (mol_id, kmer_ids) in long_read_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                pb_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { pb_vars[(var - 1) as usize].push(mol_id) }

        }
        pb_mols.insert(mol_id, kmer_ids_sorted);
    }
    for (mol_id, kmer_ids) in long_read_het_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                pb_het_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { pb_het_vars[(var - 1) as usize].push(mol_id) }

        }
        pb_het_mols.insert(mol_id, kmer_ids_sorted);
    }
    for i in 0..max_var {
        txg_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
        hic_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
        pb_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
        pb_het_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
    }
    let vars = Variants{
        paired_het_variants: paired_het_variants,
        linked_read_variants: txg_vars,
        hic_variants: hic_vars,
        long_read_variants: pb_vars,
        long_read_het_variants: pb_het_vars,
        hom_linked_read_variants: hom_linked_read_kmers,
        hom_long_read_variants: hom_long_read_kmers,
        hom_hic_variants: hom_hic_kmers,
        het_linked_read_variants: het_linked_read_kmers,
        het_long_read_variants: het_long_read_kmers,
        het_hic_variants: het_hic_kmers,
    };
    let mols = Molecules{
        longread_molid_offsets: longread_mol_id_starts,
        linked_read_molecules: txg_mols,
        hic_molecules: hic_mols,
        long_read_molecules: pb_mols,
        long_read_molecule_list: long_read_molecule_list,
        long_read_molecule_positions: long_read_molecule_positions,
        long_read_het_molecules: pb_het_mols,
        long_read_het_molecule_list: long_read_het_molecule_list,
        hom_linked_read_molecules: hom_linked_read_molecules,
        hom_long_read_molecules: hom_long_read_molecules,
        hom_hic_molecules: hom_hic_molecules,
        het_linked_read_molecules: het_linked_read_molecules,
        het_long_read_molecules: het_long_read_molecules,
        het_hic_molecules: het_hic_molecules,
    };
    (vars, mols)
}

#[derive(PartialEq)]
pub enum KmerType {
    PairedHet,
    UnpairedHet,
    Homozygous,
}

pub struct Kmers {
    kmers: HashMap<i32, String>,
    kmer_counts: HashMap<i32, i32>,
    kmer_type: HashMap<i32, KmerType>,
}

impl Kmers {
    pub fn pair(v: i32) -> i32 {
        if v < 0 { 
            if v % 2 == 0 { v + 1 } 
            else { v - 1 } 
        }
        else { 
            if v % 2 == 0 { v - 1 } 
            else { v + 1 } 
        }
    }

    pub fn load_kmers(kmerfile: &String) -> Kmers {
        let mut kmers: HashMap<i32, String> = HashMap::new();
        let mut kmer_id: i32 = 1; // no 0 kmer id as we are using sign for which of the pair
        let mut kmer_type: HashMap<i32, KmerType> = HashMap::new();
        let mut kmer_counts: HashMap<i32, i32> = HashMap::new();
        let reader = File::open(kmerfile).expect("cannot open kmer file");
        let mut reader = BufReader::new(reader);
        let mut buf1 = vec![]; let mut buf2 = vec![];
        let mut buf3 = vec![]; let mut buf4 = vec![];
        loop {
            buf1.clear(); buf2.clear(); buf3.clear(); buf4.clear();

            let bytes1 = reader.read_until(b'\t', &mut buf1).expect("cannot read file");
            if bytes1 == 0 { break; } 
    
            let bytes2 = reader.read_until(b'\t', &mut buf2).expect("cannot read file");
            if bytes2 == 0 { break; } 

            let bytes3 = reader.read_until(b'\t', &mut buf3).expect("cannot read file");
            if bytes3 == 0 { break; } 

            let bytes4 = reader.read_until(b'\n', &mut buf4).expect("cannot read file");
            if bytes4 == 0 { break; } 

            let switchme = std::str::from_utf8(&buf3).unwrap();
            match switchme {
                "HET\t" => {
                    kmer_type.insert(kmer_id, KmerType::UnpairedHet);
                    kmer_type.insert(-kmer_id, KmerType::UnpairedHet);
                    kmers.insert(kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                    //println!("")
                },
                "HOM\t" => {
                    kmer_type.insert(kmer_id, KmerType::Homozygous);
                    kmer_type.insert(-kmer_id, KmerType::Homozygous);
                },
                _ => {
                    kmer_type.insert(kmer_id, KmerType::PairedHet);
                    kmer_type.insert(-kmer_id, KmerType::PairedHet);
                    kmer_type.insert(kmer_id+1, KmerType::PairedHet);
                    kmer_type.insert(-(kmer_id+1), KmerType::PairedHet);
                    kmers.insert(kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                    kmers.insert(-kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                    let count = std::str::from_utf8(&buf2[0..(bytes2-1)]).unwrap().to_string().parse::<i32>().unwrap();
                    kmer_counts.insert(kmer_id, count);
                    kmer_counts.insert(-kmer_id, count);
                    kmers.insert(kmer_id+1, std::str::from_utf8(&buf3[0..(bytes3-1)]).unwrap().to_string());
                    kmers.insert(-(kmer_id+1), std::str::from_utf8(&buf3[0..(bytes3-1)]).unwrap().to_string());
                    let count = std::str::from_utf8(&buf4[0..(bytes4-1)]).unwrap().to_string().parse::<i32>().unwrap();
                    kmer_counts.insert(kmer_id+1, count);
                    kmer_counts.insert(-(kmer_id+1), count);
                },
            }
            kmer_id += 2;
        }
        //(kmers, kmer_counts)
        Kmers{
            kmers: kmers,
            kmer_counts: kmer_counts,
            kmer_type: kmer_type,
        }
    }
}

fn get_reader(filename: String) -> BufReader<Box<dyn Read>> {
    let filetype: Vec<&str> = filename.split(".").collect();
    let filetype = filetype[filetype.len()-1];
    let file = match File::open(filename.clone()) {
        Ok(file) => file,
        Err(error) => panic!("There was a problem opening the file: {:?}", error),
    };
    let reader: Box<dyn Read> = match filetype {
        "gz" => Box::new(GzDecoder::new(file)),
        _ => Box::new(file),
    };
    BufReader::new(reader)
}



#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
