use axgeom::*;
use axgeom::Axis;
use std::collections::BTreeMap;
use axgeom::ordered_float::*;
use duckduckgeo::grid;

use bothash::BotCollisionHash;
    mod bothash{
        use core::cmp::Ordering;
        use core::marker::PhantomData;
        #[derive(Copy,Clone)]
        pub struct BotCollisionHash<T>((usize,usize),PhantomData<T>);
        
         impl<T> PartialEq for BotCollisionHash<T> {
            fn eq(&self, other: &Self) -> bool {
                self.0.eq(&other.0)
            }
        }
        impl<T> Eq for BotCollisionHash<T> {}
        impl<T> PartialOrd for BotCollisionHash<T>{
            fn partial_cmp(&self, other: &Self) -> Option<Ordering>{
                self.0.partial_cmp(&other.0)
            }
        }
        impl<T> Ord for BotCollisionHash<T>{
            fn cmp(&self, other: &Self) -> Ordering{
                self.0.cmp(&other.0)
            }
        }
        impl<T> BotCollisionHash<T>{
            pub fn new(a:&T,b:&T)->BotCollisionHash<T>{                
                let a=a as *const _ as usize;
                let b=b as *const _ as usize;
                let [a,b]=if a<b{
                    [a,b]
                }else{
                    [b,a]
                };
                BotCollisionHash((a,b),PhantomData)
            }
        }
    }


    /*
#[derive(PartialOrd,PartialEq,Eq,Ord,Copy,Clone)]
struct BotCollisionHash(usize,usize);
impl BotCollisionHash{
    fn new<T>(a:&T,b:&T)->BotCollisionHash{                
        let a=a as *const _ as usize;
        let b=b as *const _ as usize;
        let [a,b]=if a<b{
            [a,b]
        }else{
            [b,a]
        };
        BotCollisionHash(a,b)
    }
}
*/

#[derive(PartialOrd,PartialEq,Eq,Ord,Copy,Clone)]
struct WallCollisionHash{
    a:usize,
    dir:grid::CardDir
}

fn single_hash<T>(a:&T,dir:grid::CardDir)->WallCollisionHash{
    WallCollisionHash{a:a as *const _ as usize,dir}
}



pub struct CollisionVelocitySolver<T>{
    last_bot_col:BTreeMap<BotCollisionHash<T>,f32>,
    last_wall_col:BTreeMap<WallCollisionHash,f32>
}

impl<T:Send+Sync> CollisionVelocitySolver<T>{
    pub fn new()->CollisionVelocitySolver<T>{
        CollisionVelocitySolver{last_bot_col:BTreeMap::new(),last_wall_col:BTreeMap::new()}
    }



    pub fn solve<A:Axis>(
        &mut self,
        radius:f32,
        grid_viewport:&grid::GridViewPort,
        walls:&grid::Grid2D,
        tree:&mut dinotree_alg::collectable::CollectableDinoTree<A,NotNan<f32>,T>,
        pos_func:impl Fn(&T)->&Vec2<f32> + Send + Sync +Copy,
        vel_func:impl Fn(&mut T)->&mut Vec2<f32> + Send + Sync +Copy){
        
        let diameter=radius*2.0;
        let diameter2=diameter*diameter;
        let bias_factor=0.3;
        let allowed_penetration=radius*0.2;
        let num_iterations=5;
    
        let mut collision_list={
            let ka3 = &self.last_bot_col;
            tree.collect_intersections_list_par(|a,b|{
                
                let offset=*pos_func(b)-*pos_func(a);
                let distance2=offset.magnitude2();
                if distance2>0.00001 && distance2<diameter2{
                    let distance=distance2.sqrt();
                    let offset_normal=offset/distance;
                    let separation=(diameter-distance)/2.0;
                    let bias=-bias_factor*(1.0/num_iterations as f32)*( (-separation+allowed_penetration).min(0.0));
                    
                    let hash=BotCollisionHash::new(a,b);
                    let impulse=if let Some(&impulse)=ka3.get(&hash){ //TODO inefficient to check if its none every time
                        let k=offset_normal*impulse;
                        *vel_func(a)-=k;
                        *vel_func(b)+=k;
                        impulse
                    }else{
                        0.0
                    };

                    Some((offset_normal,bias,impulse))
                }else{
                    None
                }
            })
        };

        //Package in one struct
        //so that there is no chance of mutating it twice
        #[derive(Debug)]
        struct WallCollision{
            collisions:[Option<(f32,Vec2<f32>,grid::CardDir,f32)>;2],
        }

        let mut wall_collisions={
            let ka3 = &self.last_wall_col;

            tree.collect_all_par(|rect,a|{
                
                let arr=duckduckgeo::grid::collide::is_colliding(&walls,&grid_viewport,rect.as_ref(),radius);
                let create_collision=|bot:&mut T,dir:grid::CardDir,seperation:f32,offset_normal:Vec2<f32>|{
                    let bias=-bias_factor*(1.0/num_iterations as f32)*( (-seperation+allowed_penetration).min(0.0));

                    let impulse=if let Some(&impulse)=ka3.get(&single_hash(bot,dir)){ //TODO inefficient to check if its none every time
                        let k=offset_normal*impulse;
                        *vel_func(bot)+=k;
                        impulse
                    }else{
                        0.0
                    };
                    (bias,offset_normal,dir,impulse)
                };
                match arr[0]{
                    Some((seperation,dir,offset_normal))=>{
                        
                        let wall=match arr[1]{
                            Some((seperation,dir,offset_normal))=>{
                                let seperation=seperation*2.0f32.sqrt(); //Since we are pushing diagonally dont want to over push.
                                let first=Some(create_collision(a,dir,seperation,offset_normal));
                                let second=Some(create_collision(a,dir,seperation,offset_normal));
                                WallCollision{collisions:[first,second]}
                            },
                            None=>{
                                let first=Some(create_collision(a,dir,seperation,offset_normal));
                                WallCollision{collisions:[first,None]}
                            }
                        };
                        Some(wall)
                    },
                    None=>{
                        None
                    }
                }
            })
        };

        for _ in 0..num_iterations{

            collision_list.for_every_pair_mut_par(tree,|a,b,&mut (offset_normal,bias,ref mut acc)|{
                
                let vel=*vel_func(b)-*vel_func(a);

                let mass=0.2;
                let impulse=(bias-vel.dot(offset_normal))*mass;
                
                let p0=*acc;
                *acc=(p0+impulse).max(0.0);
                let impulse=*acc-p0;
                
                let k=offset_normal*impulse;
                *vel_func(a)-=k;
                *vel_func(b)+=k;
            });     

            wall_collisions.for_every_mut_par(tree,|bot,wall|{
                
                
                for k in wall.collisions.iter_mut(){
                    if let &mut Some((bias,offset_normal,_dir,ref mut acc))=k{
                        
                        let impulse=bias-vel_func(bot).dot(offset_normal);

                        let p0=*acc;
                        *acc=(p0+impulse).max(0.0);
                        let impulse=*acc-p0;

                        *vel_func(bot)+=offset_normal*impulse;
                    }
                }; 
            })
        }  


        self.last_bot_col.clear();
        self.last_wall_col.clear();

        let (ka2,ka3):(BTreeMap<_,_>,BTreeMap<_,_>)=rayon::join(||{
            collision_list.get(&tree).iter().flat_map(|a|a.iter()).map(|(a,b,(_,_,impulse))|{
                (BotCollisionHash::new(*a,*b),*impulse)
            }).collect()
        },
        ||{
            wall_collisions.get(&tree).iter().flat_map(|(bot,wall)|{
                let k=wall.collisions.iter().filter(|a|a.is_some()).map(|a|a.unwrap());
                k.map(move |(_,_,dir,impulse)|{
                    (single_hash(bot,dir),impulse)
                })
            }).collect()
        });

        self.last_bot_col=ka2;
        self.last_wall_col=ka3;
    }
}
